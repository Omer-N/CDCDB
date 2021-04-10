import json
from tqdm import tqdm
import pandas as pd
import math

tqdm.pandas()


class ClinicalTrialsSchemaTransformer(object):
    """
    Used to create normalized and unormalized version of the df
    """

    def __init__(self):
        pass

    def transform_normalized(self, df: pd.DataFrame, dbid_to_compound_size_df: pd.DataFrame) -> dict:
        """
        :param df: raw df
        :param dbid_to_compound_size_df:
        :rtype: dict[str -> pd.DataFrame]
        :return: dictionary of normalized tables (name to df)
        """
        design_group_df = self.extract_design_groups_df(df, dbid_to_compound_size_df)
        nct_ids_of_combs_researches_df = design_group_df[['nct_id']]
        df = df.merge(nct_ids_of_combs_researches_df.drop_duplicates())
        conditions_df = self.extract_conditions_df(df)
        mesh_terms_df = self.extract_mesh_terms_df(df)
        references_df = self.extract_references_df(df)
        trials_df = self.extract_trials_df(df)
        return {'conditions_df': conditions_df,
                'mesh_terms_df': mesh_terms_df,
                'references_df': references_df,
                'trials_df': trials_df,
                'design_group_df': design_group_df}

    def extract_conditions_df(self, df: pd.DataFrame) -> pd.DataFrame:
        relevant_cols_df = df[['nct_id', 'condition_names']].drop_duplicates()
        relevant_cols_df['condition_names'] = relevant_cols_df['condition_names'].apply(
            lambda x: [] if str(x) == 'nan' else eval(x))
        res = []
        for idx, row in tqdm(relevant_cols_df.iterrows()):
            mesh_terms = row['condition_names']
            nct_id = row['nct_id']
            if mesh_terms:
                for mesh_term in mesh_terms:
                    res.append({'nct_id': nct_id, 'condition': mesh_term, 'condition_downcase': mesh_term.lower()})
        return pd.DataFrame(res).drop_duplicates()

    def extract_mesh_terms_df(self, df: pd.DataFrame) -> pd.DataFrame:
        relevant_cols_df = df[['nct_id', 'mesh_terms']].drop_duplicates()
        relevant_cols_df['mesh_terms'] = relevant_cols_df['mesh_terms'].apply(
            lambda x: [] if str(x) == 'nan' else eval(x))
        res = []
        for idx, row in tqdm(relevant_cols_df.iterrows()):
            mesh_terms = row['mesh_terms']
            nct_id = row['nct_id']
            if mesh_terms:
                for mesh_term in mesh_terms:
                    res.append({'nct_id': nct_id, 'mesh_term': mesh_term, 'mesh_terms_downcase': mesh_term.lower()})
        return pd.DataFrame(res).drop_duplicates()

    def extract_references_df(self, df: pd.DataFrame) -> pd.DataFrame:
        relevant_cols_df = df[['nct_id', 'refs']].drop_duplicates()
        relevant_cols_df['refs'] = relevant_cols_df['refs'].apply(lambda x: [] if str(x) == 'nan' else eval(x))
        res = []
        for idx, row in tqdm(relevant_cols_df.iterrows()):
            references = row['refs']
            nct_id = row['nct_id']
            if references:
                for reference in references:
                    res.append({'nct_id': nct_id, 'reference_type': reference[0], 'reference': reference[1]})
        return pd.DataFrame(res)

    def extract_trials_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[['nct_id', 'study_start_date', 'overall_status', 'phase', 'completion_date',
                   'enrollment', 'enrollment_type', 'number_of_arms', 'number_of_groups',
                   'why_stopped']].drop_duplicates()

    def eval_and_json_double(self, x):
        return json.dumps([eval(y) for y in eval(x)])

    def eval_and_json(self, x):
        return json.dumps(eval(x))

    def extract_design_groups_df(self, df, dbid_to_compound_size_df):
        df['drugbank_identifier'] = df['identifiers_entity'].apply(lambda x: eval(x)[0])
        df['pubchem_identifier'] = df['identifiers_entity'].apply(lambda x: eval(x)[1])
        df = df.merge(dbid_to_compound_size_df, left_on="drugbank_identifier", right_on="id", how='left')
        df['is_complex_compound'] = df['compound_size'].apply(lambda x: math.isnan(x) or x > 2)
        drugbank_nutraceuticals_df = pd.read_excel("input_data/drugbank_nutraceuticals.xlsx")
        df = df.merge(drugbank_nutraceuticals_df[['Nutraceutical', 'DrugBank ID']], left_on="drugbank_identifier",
                      right_on="DrugBank ID", how='left')
        df['notNutraceutical'] = df['Nutraceutical'].apply(lambda x: math.isnan(x) or x == False)

        intervention_names_for_group = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['interventions_names'].progress_apply(list))

        selected_name_group = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['selected_name'].progress_apply(list))

        drugbank_identifiers_in_group = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['drugbank_identifier'].progress_apply(list))

        pubchem_identifiers_in_group = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['pubchem_identifier'].progress_apply(list))

        is_valid_compound_in_group = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['is_complex_compound'].progress_apply(list))

        is_contain_nutraceuticals = pd.DataFrame(
            df.groupby(['nct_id', 'design_group_id'])['notNutraceutical'].progress_apply(list))

        result_df = df[['nct_id', 'design_group_id', 'group_type', 'title']].merge(intervention_names_for_group,
                                                                                   on='design_group_id') \
            .merge(is_valid_compound_in_group, on='design_group_id') \
            .merge(is_contain_nutraceuticals, on='design_group_id') \
            .merge(selected_name_group, on='design_group_id').merge(drugbank_identifiers_in_group,
                                                                    on='design_group_id').merge(
            pubchem_identifiers_in_group, on='design_group_id').astype(str).drop_duplicates()
        result_df['drugbank_identifier'] = result_df['drugbank_identifier'].apply(self.eval_and_json)
        result_df['pubchem_identifier'] = result_df['pubchem_identifier'].apply(self.eval_and_json)
        result_df['interventions_names'] = result_df['interventions_names'].apply(self.eval_and_json_double)
        result_df['selected_name'] = result_df['selected_name'].apply(
            self.eval_and_json_double)
        result_df = result_df[result_df['is_complex_compound'].apply(lambda x: all(eval(x)))]
        result_df = result_df[result_df['notNutraceutical'].apply(lambda x: all(eval(x)))]
        result_df = result_df.drop(['is_complex_compound', 'notNutraceutical'], axis=1)

        def is_comb_contain_2_non_placebo(ids_arr):
            is_placebo = list(map(lambda y: y == 'PLACEBO', eval(ids_arr)))
            if len(is_placebo) - sum(is_placebo) >= 2:
                return True
            return False

        result_df = result_df[result_df["drugbank_identifier"].apply(is_comb_contain_2_non_placebo)]
        result_df = result_df[result_df['drugbank_identifier'].apply(
            lambda x: True if len(eval(x)) > 2 else len(eval(x)) == len(set(eval(x))))]
        return result_df


if __name__ == '__main__':
    import argparse

    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("input_path", help="path_to_raw_df.csv")
    argument_parser.add_argument("output_dir", help="output directory for the transformed tables")
    argument_parser.add_argument("dbid_to_compound_size_df", default="input_data/dbid_to_compound_size_df.csv",
                                 help="dbid to compound size df path")
    args = argument_parser.parse_args()
    raw_df = pd.read_csv(args.input_path)
    dbid_to_compound_size_df = pd.read_csv(args.dbid_to_compound_size_df)
    clinical_trials_schema_transformer = ClinicalTrialsSchemaTransformer()
    normalized_tables = clinical_trials_schema_transformer.transform_normalized(raw_df, dbid_to_compound_size_df)
    for name, df in normalized_tables.items():
        df.to_csv(f'{args.output_dir}/{name}.csv', index=False)
