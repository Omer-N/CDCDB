# CDCDB
To create a new version
1. create python venv
2. fill `drug_combs/input_data/aact_credentials.json` with your credentials for AACT
3. > pip install -r requirements
4. > cd drug_combs
5. Create a new version by running the creation script 
  > ./create_version.sh
6. Follow the progress-bars, and install scispacy's model if missing


Note, there are plenty of caches used in this system, therefore the first run would be longer.
