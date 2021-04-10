import json
from sqlalchemy import create_engine
import pandas as pd
import logging
import sys

'''
Number of studies by year from AACT
SELECT EXTRACT(YEAR from start_date), COUNT(nct_id) FROM studies GROUP BY EXTRACT(YEAR from start_date)
'''


class AACTFetcher(object):
    """
    This class is used to retrive latest clinical-trials info from AACT (which is updated every 24 hours)
    """

    def __init__(self, aact_url, aact_db_username, aact_db_password):
        self.db_connection = self.get_aact_db_connection(aact_url, aact_db_username, aact_db_password)

    def get_aact_db_connection(self, aact_url, aact_db_username, aact_db_password):
        logging.info("Initiating db connection")
        engine = create_engine(
            f'postgresql://{aact_db_username}:{aact_db_password}@{aact_url}/aact')
        connection = engine.connect()
        return connection

    def fetch_data_frame(self):
        logging.info("Fetching dataframe from remote")
        return pd.read_sql(self.get_query(), self.db_connection)

    def get_query(self):
        logging.log(logging.DEBUG, "query requested")

        return '''
        SELECT relevant_studies.*,
       collected_conditions.mesh_terms,
       collected_conditions.downcase_mesh_terms,
       collected_conditions.condition_names,
       collected_conditions.condition_downcase_names
FROM (SELECT studies.nct_id                                                       nct_id,
             studies.start_date                                                   study_start_date,
             studies.completion_date                                              completion_date,
             studies.enrollment                                                   enrollment,
             studies.enrollment_type                                              enrollment_type,
             studies.number_of_arms                                               number_of_arms,
             studies.number_of_groups                                             number_of_groups,
             studies.why_stopped                                                  why_stopped,
             studies.phase                                                        phase,
             studies.overall_status                                               overall_status,
             studies.last_known_status                                            last_known_status,
             studies.is_fda_regulated_drug                                        is_fda_regulated_drug,
             dg.id                                                                design_group_id,
             dgi.intervention_id                                                  interventions_id,
             dg.group_type,
             dg.title,
             i.name                                                               intervention_names,
             CASE
                 WHEN ions.intervention_other_names is null THEN json_build_array(i.name, json_build_array())
                 ELSE json_build_array(i.name, ions.intervention_other_names) END interventions_with_other_names,
             i.description                                                        intervention_description,
             collected_refs.refs                                                  refs

      FROM studies
               LEFT JOIN design_groups dg
                         on studies.nct_id = dg.nct_id
               LEFT JOIN design_group_interventions dgi on dg.id = dgi.design_group_id
               LEFT JOIN interventions i on dgi.intervention_id = i.id
               LEFT JOIN (SELECT json_agg(ion.name) intervention_other_names, intervention_id
                          FROM intervention_other_names ion
                          GROUP BY intervention_id) ions on i.id = ions.intervention_id
               LEFT JOIN (SELECT nct_id,
                                 json_agg(json_build_array(study_references.reference_type,
                                                           study_references.citation)) refs
                          FROM study_references
                          GROUP BY nct_id
      ) as collected_refs on collected_refs.nct_id = studies.nct_id
      where intervention_type = 'Drug'
        and dg.id in (SELECT dg.id
                      FROM studies
                               LEFT JOIN design_groups dg on studies.nct_id = dg.nct_id
                               LEFT JOIN design_group_interventions dgi on dg.id = dgi.design_group_id
                               LEFT JOIN interventions i on dgi.intervention_id = i.id
                      where intervention_type = 'Drug'
                      GROUP BY studies.nct_id
                             , dg.id
                      HAVING COUNT(*)
                                 > 1)
     ) relevant_studies
         LEFT OUTER JOIN (
    SELECT collected_conds.nct_id nct_id,
           mesh_terms,
           downcase_mesh_terms,
           condition_names,
           condition_downcase_names
    FROM (SELECT bc.nct_id                       nct_id,
                 json_agg(bc.mesh_term)          mesh_terms,
                 json_agg(bc.downcase_mesh_term) downcase_mesh_terms
          FROM browse_conditions bc
          GROUP BY bc.nct_id) bc_collected
             LEFT JOIN (SELECT c.nct_id                  nct_id,
                               json_agg(c.name)          condition_names,
                               json_agg(c.downcase_name) condition_downcase_names
                        FROM conditions c
                        GROUP BY c.nct_id) collected_conds
                       ON collected_conds.nct_id = bc_collected.nct_id
) collected_conditions ON relevant_studies.nct_id = collected_conditions.nct_id;
'''
    def close_connection(self):
        self.db_connection.close()

def read_creds_from_file(path):
    try:
        credentials_file = open(path)
        cred = json.load(credentials_file)
        return cred
    except IOError as e:
        print(f"Failed to open AACT credentials file: {e}")
        raise e


if __name__ == '__main__':
    argv = sys.argv
    if len(argv) != 3:
        sys.exit("Usage: aact_fetcher.py <aact_credentials_file.json> <output_path.csv>")

    creds = read_creds_from_file(argv[1])
    aact_fetcher = AACTFetcher(creds['url'], creds['username'], creds['password'])
    df = aact_fetcher.fetch_data_frame()
    print(f'Fetched total of: {len(df)} rows')
    output_path = argv[2]
    df.to_csv(output_path, index=False)
    print(f'Saved successfully the fetched dataframe to {output_path}')
    # pip install psycopg2-binary
