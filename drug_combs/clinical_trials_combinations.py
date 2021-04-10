import sys

sys.path.insert(0, '../..')

import edlib
from functools import lru_cache
from src.drug_combs.aact_fetcher import AACTFetcher
from src.drug_combs.add_identifier_to_df import DataframeDrugIdentifiersAdder
from src.drug_identfiers_resolver.identifiers_resolver import *
import re
import json
import argparse
import warnings
import diskcache as dc

import spacy
from scispacy.linking import EntityLinker

nlp = spacy.load("en_core_sci_lg")
linker = EntityLinker(resolve_abbreviations=True, name="umls")

nlp.add_pipe(linker)

warnings.filterwarnings("ignore")
DRUG_IDENTIFIERS_COLUMN = "drug_identifiers"

INTERVENTIONS_WITH_OTHER_NAMES_COL = 'interventions_with_other_names'

INTERVENTIONS_NAMES_COL = 'interventions_names'
INTERVENTIONS_NAMES_CLEANED_COL = 'interventions_names_cleaned'

drugs_TUIs = {
    "T109", "T114", "T116", "T121", "T123", "T125", "T126", "T129", "T195", "T200"
}
# not_allowed_TUIs = {"T028", "T073", "T061"}
# not_allowed_TUIs = {"T127", "T197"}


class DataPreProcessor(object):
    """
    Used to clean the data and add relevant fields like drug identifier
    """

    def __init__(self):
        pass

    def _process(self, raw_df, pipeline_functions):
        result_df = raw_df
        for processing_function in pipeline_functions:
            result_df = processing_function(result_df)
        return result_df


class AACTDataPreProcessor(DataPreProcessor):

    def __init__(self, drug_identifiers_resolver: DrugIdentifiersResolver, drug_resolving_threads=16,
                 ner_cache_path='NER-mappings.cache'):
        super().__init__()
        self.drug_identifiers_resolver = drug_identifiers_resolver
        self.cache = dc.Cache(ner_cache_path)

    def preprocess(self, raw_df: pd.DataFrame) -> pd.DataFrame:
        functions_pipe = [
            self.flatten_interventions,
            self.clean_drug_names,
            self.clean_placebos,
            self.extract_entities,
            # self.add_identifiers, the problem is with multithreads including the language model
            # self.split_identifiers
        ]
        return self._process(raw_df, functions_pipe)

    def remove_special_words(self, x):
        words_to_remove = ["single", "dose", "low", "slow", "soft tissue", "solid", "spray", "stable", "subarachnoid",
                           "subconjunctival", "subcutaneous", "sublingual", "sublingual", "submucosal", "suppositories",
                           "sustained-release", "tablet", "tablets", "therapy", "topical", "transdermal",
                           "transmucosal", "transplacental", "transtracheal", "transtympanic", "treatment", "troches",
                           "ureteral", "urethral", "usual care", "vaginal", "%", "mg", "kg", "mg/day", "oral",
                           "suspension", "low dose", "fixed", "combination", "drops"]
        words = x.split()
        result = ' '.join([word for word in words if word.lower() not in words_to_remove])
        return result.strip()

    def clean_drug_name(self, x, regexes):
        x = str(x)
        x = self.remove_special_words(x)
        x = re.sub('"|\'|mg/day|,|•|™|®|(?i)oral\\b|(?i)IV\\b', '', x).strip()
        x = re.sub('α', 'alfa', x).strip()
        x = re.sub('[μμμµ]', 'μ', x).strip()
        x = re.sub('[-+]?\d*\.?\d*%', '', x).strip()
        x = re.sub('([-+]?\d*\.\d*)', '', x).strip()
        for regex in regexes:
            res = regex.findall(x)
            if res != []:
                x = res[0]
        return x.strip()

    def clean_drug_names(self, df: pd.DataFrame) -> pd.DataFrame:
        comparator_regex = re.compile('Comparator: (.*)')
        remove_mg_kg = re.compile('^(.*?)(?:(?:\/\d)|(?: \d)|(?:,(?:.*)\d)).*?(?:mg|kg|μg|mcg)(?:.*?)$')
        remove_percentage_from_start = re.compile('[-+]?\d*\.?\d*%(?: )?(.*?)$')
        remove_dosage_before = re.compile('^\d.*(?:mg|kg|μg|mcg|µg)(.*?)$')
        remove_par = re.compile('^(.*?)\(.*?\).*$')
        remove_rect_par = re.compile('^(.*?)\[.*?\].*$')
        regs = [comparator_regex, remove_mg_kg, remove_percentage_from_start, remove_dosage_before, remove_par,
                remove_rect_par]
        # logging.info("cleaning drug names")
        df[INTERVENTIONS_NAMES_CLEANED_COL] = df[INTERVENTIONS_NAMES_COL].progress_apply(
            lambda names: [self.clean_drug_name(name, regs) for name in names])
        return df

    def flatten_interventions(self, df: pd.DataFrame) -> pd.DataFrame:
        # logging.info("flattening interventions array")
        df[INTERVENTIONS_NAMES_COL] = df[INTERVENTIONS_WITH_OTHER_NAMES_COL].astype(str).apply(self.flatten_array)
        return df

    def replace_placebo(self, interventions_names):
        if 'placebo' in str(interventions_names).lower():
            return ['placebo']
        return interventions_names

    def clean_placebos(self, df: pd.DataFrame) -> pd.DataFrame:
        # logging.info("flattening interventions array")
        df[INTERVENTIONS_NAMES_CLEANED_COL] = df[INTERVENTIONS_NAMES_CLEANED_COL].apply(self.replace_placebo)
        return df

    def add_identifiers(self, df: pd.DataFrame) -> pd.DataFrame:
        adder = DataframeDrugIdentifiersAdder(self.drug_identifiers_resolver)
        return adder.add_identifiers_column(df, INTERVENTIONS_NAMES_CLEANED_COL, DRUG_IDENTIFIERS_COLUMN, 16)

    def split_identifiers(self, df: pd.DataFrame) -> pd.DataFrame:
        df['drugbank_identifier'] = df[DRUG_IDENTIFIERS_COLUMN].apply(lambda x: x[0])
        df['pubchem_identifier'] = df[DRUG_IDENTIFIERS_COLUMN].apply(lambda x: x[1])
        return df

    def flatten_array(self, arr):
        arr = eval(arr)
        return [arr[0]] + arr[1]

    def extract_entities(self, df: pd.DataFrame):
        selected_name = df[INTERVENTIONS_NAMES_CLEANED_COL].progress_apply(self.extract_entities_from_list)
        df['selected_name'] = selected_name
        errors = df[df['selected_name'].apply(
            lambda x: x == [] or x == 'ERROR: contained more than one entity-should drop')]
        errors.to_csv("errors.csv")
        df = df[~df['design_group_id'].isin(errors['design_group_id'].unique())]
        return df

    @lru_cache(50000)
    def get_most_relevant_name(self, ent):
        min_alias_val = 9999999
        min_alias = None
        ent_str = str(ent)
        all_types = []
        for x in ent._.kb_ents:
            entity = linker.kb.cui_to_entity[x[0]]
            all_types += entity.types
            contains_at_least_one_require_type = len(drugs_TUIs - set(entity.types)) != len(drugs_TUIs)
            if contains_at_least_one_require_type:
                for alias in entity.aliases:
                    distance = edlib.align(alias, ent_str)['editDistance']
                    if distance < min_alias_val:
                        min_alias_val = distance
                        min_alias = entity.canonical_name

        # contains_unwanted_tui = len(not_allowed_TUIs) != len(not_allowed_TUIs - set(all_types))
        # if contains_unwanted_tui:
        #     return None
        return min_alias

    def resolve_details(self, list_of_drugs):
        res = []
        for drug_syns in list_of_drugs:
            syns_res = []
            for syn in drug_syns:
                if syn:
                    umls_concepts = syn[0].ents[0]._.umls_ents
                    for concept in umls_concepts:
                        resolved_concept = linker.kb.cui_to_entity[concept[0]]
                        syns_res.append(resolved_concept)
            res.append(syns_res)
        return res

    def extract_entities_from_list(self, drugs_list):
        if drugs_list == ['placebo']:
            return drugs_list
        res = []
        for x in drugs_list:
            name = self.cache.get(x)
            if name is not None:
                res.append(name)
                continue
            doc = nlp(str(x))
            if doc.ents:
                if len(doc.ents) > 1:
                    return "ERROR: contained more than one entity-should drop"
                name = self.get_most_relevant_name(doc.ents[0])
                if name is not None:
                    self.cache[x] = name
                    res.append(name)
        return res


class DatasetCreator(object):
    """
    Used to generate the dataset of all the combinations in clinical trials, from file or from AACT DB.
    """

    def __init__(self, qid_to_drugbank_path='input_data/qid_to_drugbank.json',
                 qid_to_pubchem_path='input_data/qid_to_pubchem.json', wikidata_cache_path='wikidata_cache.csv',
                 api_cache_path=''):
        self.qid_to_drugbank_path = qid_to_drugbank_path
        self.qid_to_pubchem_path = qid_to_pubchem_path
        self.wiki_cache_path = wikidata_cache_path
        self.api_cache_path = api_cache_path

    def process_df(self, df):
        wikidata_ids_resolver = WikiDataIdsResolver(self.qid_to_drugbank_path, self.qid_to_pubchem_path)
        resolver = DrugIdentifiersResolver(wikidata_ids_resolver, APIBasedIdentifiersResolver(self.api_cache_path))
        aact_preprocessor = AACTDataPreProcessor(resolver)
        return aact_preprocessor.preprocess(df)

    def create_updated_dataset(self, aact_url, aact_db_username, aact_db_password):
        aact_fetcher = AACTFetcher(aact_url, aact_db_username, aact_db_password)
        up_to_date_df = aact_fetcher.fetch_data_frame()
        aact_fetcher.close_connection()
        return self.process_df(up_to_date_df)


def main(args):
    dataset_creator = DatasetCreator()
    if args.input_path is not None:
        input_file = pd.read_csv(args.input_path)[:100]
        processed_df = dataset_creator.process_df(input_file)
    else:
        try:
            path = args.aact_params_file_path
            if path is None:
                exit("You must provide aact credentials file or input file")
            credentials_file = open(path)
            cred = json.load(credentials_file)
            processed_df = dataset_creator.create_updated_dataset(cred['url'], cred['username'], cred['password'])
        except IOError as e:
            print(f"Failed to open AACT credentials file: {e}")
    processed_df.to_csv(args.output_path, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--wiki-cache", help="path to the wikidata api cache json file")
    parser.add_argument("--api-cache", help="path to the drugbank/pubchem api cache csv file")
    parser.add_argument("--qid-to-dbid", help="path to the qid to dbid mapping file")
    parser.add_argument("--qid-to-pubchem", help="path to the qid to pubchem_id mapping file")
    parser.add_argument("--quiet", default=False)
    parser.add_argument("--input_path", default=None, type=str,
                        help="AACT raw dataframe path, if not provided then up to date snapshot is fetched")
    parser.add_argument("--aact_params_file_path", type=str,
                        help="path to file contains the db url, name and password of AACT")
    parser.add_argument("output_path", type=str, help="path to write the output csv")
    args = parser.parse_args()
    main(args)
    # usage example python clinical_trials_combinations.py --input_path aact_unaggregated_data.csv data/clinical_trials_comb.csv
