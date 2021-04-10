import diskcache as dc

import pubchempy as pcp
import pandas as pd
from os import path
import json
import urllib

import requests

CACHE_SEP = "@@@@@@"

MISSING_PUBCHEM_ID = '-1'
PLACEBO_CODE = "PLACEBO"
MISSING_DRUGBANK_ID = '-1'

DRUG_BANK_NAME_SEARCH_URL = f'https://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&query=drug_name&searcher=drugs'
PUBCHEM_SEARCH_BY_NAME_URL = \
    f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/drug_name/xrefs/RegistryID,RN,PubMedID/JSONP'
PUBCHEM_SEARCH_BY_CID_URL = \
    f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/cid_value/xrefs/RegistryID,RN,PubMedID/JSONP'


class APIBasedIdentifiersResolver(object):
    def __init__(self, path_to_cache=''):
        self.csvs_dir = path_to_cache
        csv_path = f'dbid_disk_cache'
        self.cache = dc.Cache(csv_path)

    def get_sids_by_name(self, drug_name):
        drug_name = self.process_query(drug_name)
        csv_path = self.csvs_dir + 'sides.csv'
        exists, result = self.check_query_in_file(csv_path, drug_name)
        if exists:
            return result
        else:
            result = pcp.get_sids(drug_name, 'name')
            sids = result[0]['SID']
            self.add_query_to_file(drug_name, sids)
            return sids

    def get_aids_by_name(self, drug_name):
        drug_name = self.process_query(drug_name)
        csv_path = self.csvs_dir + 'aids.csv'
        exists, result = self.check_query_in_file(csv_path, drug_name)
        if exists:
            return result
        else:
            result = pcp.get_aids(drug_name, 'name')
            aids = result[0]['AID']
            self.add_query_to_file(drug_name, aids)
            return aids

    def get_cids_by_name(self, drug_name):
        drug_name = self.process_query(drug_name)
        csv_path = self.csvs_dir + 'cides.csv'
        exists, result = self.check_query_in_file(csv_path, drug_name)
        if exists:
            return result
        else:
            result = pcp.get_cids(drug_name, 'name')
            cids = result[0]['CID']
            self.add_query_to_file(drug_name, cids)
            return cids

    def get_all_ids(self, drug_name):
        aids = self.get_aids_by_name(drug_name)
        print(aids)
        sids = self.get_sids_by_name(drug_name)
        print(sids)
        # cids = get_cids_by_name(drug_name)
        return aids, sids

    def get_cid_and_dbid_by_name(self, drug_name):
        """
        Return both CID and DrugBank ID of a drug
        :param drug_name: drug name
        :return: Dictionary that contains CID (key is 'CID') and DrugBank ID (key is 'DB_ID')
        """
        drug_name = self.process_query(drug_name)
        url = PUBCHEM_SEARCH_BY_NAME_URL.replace('drug_name', str(drug_name))
        try:
            response = urllib.request.urlopen(url)
            json_string = response.read().decode()[9:-3]
            data = json.loads(json_string)
            registry_ids = data['InformationList']['Information'][0]['RegistryID']
            cid = data['InformationList']['Information'][0]['CID']
            for id in registry_ids:
                if id[:2] == 'DB':
                    return {'CID': cid, 'DB_ID': id}
        except urllib.error.HTTPError as e:
            print(f"HTTP error raised during handling of {drug_name}: {e}")

    @staticmethod
    def get_drug_bank_id_by_cid(cid):
        """
        Returns the DrugBank ID of a given CID
        :param cid: CID
        :return: DrugBank ID
        """
        url = PUBCHEM_SEARCH_BY_CID_URL.replace('cid_value', str(cid))
        response = urllib.request.urlopen(url)
        json_string = response.read().decode()[9:-3]
        data = json.loads(json_string)
        registry_ids = data['InformationList']['Information'][0]['RegistryID']
        for id in registry_ids:
            if id[:2] == 'DB':
                return id

    def get_drug_bank_code_by_name(self, drug_name):
        """
        Returns drugbank drug code of a given drug.
        The function check first if the drug code exists in the saved file, if not checks with drugbank
        :param drug_name: name of the drug to search
        :return: drugbank drug code
        """
        drug_name = self.process_query(drug_name)
        csv_path = self.csvs_dir + 'drugbank_codes.csv'
        exists, result = self.check_query_in_file(csv_path, drug_name)
        # exists = False
        if exists:
            return result
        else:
            try:
                drug_code = self._retrieve_from_drugbank(drug_name)
                if drug_code[:2] == 'DB':
                    self.add_query_to_file(drug_name, drug_code)
                    return drug_code
                else:
                    self.add_query_to_file(drug_name, 'drug not found')
                    return 'drug not found'
            except urllib.error.HTTPError as e:
                print(f"encounter {e} while trying to fetch {drug_name}")
            except Exception as e:
                print(f"Encountered unhandled exception {e} while fetching {drug_name}, ignoring error returning None")

    def _retrieve_from_drugbank(self, drug_name):
        """
        Search drugbank for the drug code
        :param drug_name: drug to be searched
        :return: drug code
        """
        drug_name = self.process_query(drug_name)
        url = DRUG_BANK_NAME_SEARCH_URL.replace('drug_name', str(drug_name))
        response = requests.get(url)
        drug_code = response.url[response.url.rfind('/') + 1:]
        return drug_code

    @staticmethod
    def process_query(query):
        modified_query = query.strip()
        modified_query = modified_query.strip('\t')
        modified_query = modified_query.replace(' ', '+')
        return modified_query

    def check_query_in_file(self, file_path, query):
        dbid_in_cache = self.cache.get(query)
        if dbid_in_cache is None:
            return False, None
        return True, dbid_in_cache

    def add_query_to_file(self, query, result):
        self.cache[query] = str(result)


class WikiDataIdsResolver(object):
    def __init__(self, qid_to_dbid_path, qid_to_pubchem_id_path, cache_file_path='default_wikicache'):
        """
        This object fetches identifiers for drugbank and pubchem from wikidata given a name by utilizing Wikidata's
        search engine's API
        :param qid_to_dbid_path: a path to a CSV of mapping between QIDS (wikidata ids) to drugbank IDs
        (might be generated using the following sparql query: `SELECT * WHERE { SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". } OPTIONAL { ?item wdt:P715 ?_____Drugbank. } }` )
        :param qid_to_pubchem_id_path: a path to a CSV of mapping between QIDS (wikidata ids) to drugbank IDs
         (might be generated using the following sparql query: `SELECT * WHERE { SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". } OPTIONAL { ?item wdt:P662 ?__Pubchem. } }` )
        """
        print("creating wiki resolver")
        self.cache_file_path = cache_file_path
        # load qid to id mapping
        self.qid_to_dbid = self.get_qid_to_id_dict(qid_to_dbid_path)
        self.qid_to_pubchem_id = self.get_qid_to_id_dict(qid_to_pubchem_id_path)
        self.cache = dc.Cache(self.cache_file_path)

    def get_ids_by_name(self, name):
        """
        :param name: name of drug
        :return: (drugbank_id, pubchem_id), if found more than one return the first
        """
        exists, dbid, pbid = self.check_query_in_file(name)
        if exists:
            return dbid, pbid
        qids = self.get_qids_for(name)
        dbids = [self.qid_to_dbid.get(qid, MISSING_DRUGBANK_ID) for qid in qids]
        pubchem_ids = [self.qid_to_pubchem_id.get(qid, MISSING_PUBCHEM_ID) for qid in qids]
        dbid, pubchem_id = MISSING_DRUGBANK_ID, MISSING_PUBCHEM_ID
        if dbids:
            dbid = dbids[0]
            if dbid != MISSING_DRUGBANK_ID:
                dbid = "DB" + str(dbid)
        if pubchem_ids:
            pubchem_id = pubchem_ids[0]
            if pubchem_id != MISSING_PUBCHEM_ID:
                pubchem_id = "CID" + str(pubchem_id)
        self.add_query_to_file(name, dbid, pubchem_id)
        return dbid, pubchem_id

    def check_query_in_file(self, query):
        cached_entry = self.cache.get(query)
        if cached_entry is None:
            return False, None, None
        pbid, dbid = cached_entry.split(CACHE_SEP)
        return True, pbid, dbid

    def add_query_to_file(self, query, pubchem_id, dbid):
        self.cache[query] = f"{pubchem_id}{CACHE_SEP}{dbid}"

    def get_qids_for(self, search_query):
        """
        Search for wikidata pages for this name
        :param search_query: search_query
        :return: list of relevant QIDS
        """
        try:
            api_url = f"https://www.wikidata.org/w/api.php?action=wbsearchentities&search={search_query}&language=en&format=json&limit=50"
            response = requests.get(api_url)
            if response.status_code != 200:
                return '-1'
            query_results = response.json().get('search', [])
            qids = []
            for query_result in query_results:
                qids.append(query_result['id'])
            return qids
        except Exception as e:
            print(f"unexpected error {e} happened while searching for {search_query}")
            return []

    def get_qid_to_id_dict(self, path):
        with open(path) as dict_file:
            return json.load(dict_file)


class DrugIdentifiersResolver(object):
    def __init__(self, wikidata_ids_resolver: WikiDataIdsResolver,
                 api_identifiers_resolver: APIBasedIdentifiersResolver):
        self.api_identifiers_resolver = api_identifiers_resolver
        self.wikidata_ids_resolver = wikidata_ids_resolver

    def resolve_array(self, names) -> tuple:
        """
        :param names: possible names for the same drug
        :return: tuple of the (drugbank_id, pubchem_id)  found for the drug in the list
        """
        if names == '':
            return '', ''
        if names is None:
            return MISSING_DRUGBANK_ID, MISSING_PUBCHEM_ID
        result_drugbank_id, result_pubchem_id = MISSING_DRUGBANK_ID, MISSING_PUBCHEM_ID
        assert isinstance(names, list), "should be list"
        for name in names:
            name = str(name)
            if "placebo" in name.lower():
                return PLACEBO_CODE, PLACEBO_CODE
            drugbank_id, pubchem_id = self.wikidata_ids_resolver.get_ids_by_name(name)
            if result_drugbank_id == MISSING_DRUGBANK_ID and drugbank_id != MISSING_DRUGBANK_ID:
                result_drugbank_id = drugbank_id
            if result_pubchem_id == MISSING_PUBCHEM_ID and pubchem_id != MISSING_PUBCHEM_ID:
                result_pubchem_id = pubchem_id
            code_from_drugbank = self.api_identifiers_resolver.get_drug_bank_code_by_name(name)
            if code_from_drugbank != 'drug not found' and code_from_drugbank is not None:
                result_drugbank_id = code_from_drugbank
            if drugbank_id != MISSING_DRUGBANK_ID and pubchem_id != MISSING_PUBCHEM_ID:
                break

        return result_drugbank_id, result_pubchem_id


if __name__ == '__main__':
    # Example:
    id_converter = APIBasedIdentifiersResolver()
    print(id_converter.get_cid_and_dbid_by_name('aspirin'))
    print(id_converter.get_drug_bank_code_by_name('aspirin'))
    print(id_converter.get_sids_by_name('aspirin'))
    print(id_converter.get_drug_bank_id_by_cid(2244))
