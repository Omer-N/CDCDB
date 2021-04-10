import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
import sys
sys.path.insert(0, '../..')
from src.drug_identfiers_resolver.identifiers_resolver import WikiDataIdsResolver, DrugIdentifiersResolver, \
    APIBasedIdentifiersResolver
import warnings

warnings.filterwarnings("ignore")
tqdm.pandas()


class DataframeDrugIdentifiersAdder(object):
    def __init__(self, identifiers_resolver, array_like=True, as_str_array=False):
        self.identifiers_resolver = identifiers_resolver
        self.array_like = array_like
        self.as_str_array = as_str_array

    def _parallelize_dataframe(self, df, src_col, dest_col, n_cores):
        # logging.info(f"Parallel editing df from in_col: {src_col} to {dest_col}")
        self.source_col, self.dest_col = src_col, dest_col
        df_split = np.array_split(df, n_cores)
        pool = Pool(n_cores)
        mapper = self.add_identifiers_mapper  # self._map_split(src_col, dest_col, self.identifiers_resolver)
        print("Start mapping")
        df = pd.concat(pool.map(mapper, df_split))
        # df = mapper(df)
        pool.close()
        pool.join()
        return df

    def add_identifiers_column(self, df, source_col, dest_col, n_cores=10):
        return self._parallelize_dataframe(df, source_col,
                                           dest_col, n_cores)

    def add_identifiers_mapper(self, df):
        """
        :param df:
        :param source_col: array like column with possible names for drugs
        :param dest_col:
        :return: df with the column
        """
        if self.as_str_array:
            df[self.dest_col] = df[self.source_col].progress_apply(
                lambda x: self.identifiers_resolver.resolve_array(eval(x)))
        elif self.array_like:
            df[self.dest_col] = df[self.source_col].progress_apply(
                lambda x: self.identifiers_resolver.resolve_array(x))
        elif not self.array_like:
            df[self.dest_col] = df[self.source_col].progress_apply(
                lambda x: self.identifiers_resolver.resolve_array([x]))
        elif self.array_like == '/':
            df[self.dest_col] = df[self.source_col].progress_apply(
                lambda x: self.identifiers_resolver.resolve_array(x.split("/")))
        else:
            exit("Unsupported transformation")
        return df


def main(args):
    is_csv = args.input_path.endswith(".csv")
    is_xlsx = args.input_path.endswith(".xlsx")
    if is_csv:
        df = pd.read_csv(args.input_path)
    elif is_xlsx:
        df = pd.read_excel(args.input_path)
    else:
        exit("Unsupported file format")
    wikidata_ids_resolver = WikiDataIdsResolver("../drug_combs/input_data/qid_to_drugbank.json",
                                                "../drug_combs/input_data/qid_to_pubchem.json",
                                                cache_file_path="wikidata_disk_cache")
    resolver = DrugIdentifiersResolver(wikidata_ids_resolver, APIBasedIdentifiersResolver())
    drug_identifiers_adder = DataframeDrugIdentifiersAdder(resolver, args.aslist, args.as_str_array)
    print("Start resolving")
    try:
        df = drug_identifiers_adder.add_identifiers_column(df, args.input_column, args.result_column)#, args.processes)
        print(f"Resolved all: {len(df)} results")
        if is_csv:
            df.to_csv(args.output_path)
        else:
            df.to_excel(args.output_path)
        print(f"Saved results to {args.output_path}")
    except Exception as e:
        print(f"Failed in resolving: {e}")

if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("input_path", type=str, help="Path to input csv/xlsx file")
    argument_parser.add_argument("input_column", type=str, help="The column in the file that contains the names")
    argument_parser.add_argument("result_column", type=str,
                                 help="The column in the CSV that should contain the identifiers")
    argument_parser.add_argument("--aslist", type=bool, help="Whether the input column is list of drugs or single drug",
                                 default=False)
    argument_parser.add_argument("--as_str_array", type=bool, help="Whether the input column is list of drugs or single drug",
                                 default=False)
    argument_parser.add_argument("--processes", default=10, type=bool, help="number of processes to use")
    argument_parser.add_argument("output_path", type=str, help="The path to save the result csv")
    args = argument_parser.parse_args()
    main(args)

