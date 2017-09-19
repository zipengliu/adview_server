import os

# developement or production
ENVIRONMENT = 'development'

MONGODB_URL_DEV = 'mongodb://localhost:27017/visphy_dev'
MONGODB_URL_PROD = 'mongodb://visphybe.cs.ubc.ca:27017/visphy_dev'

# compress the json returned to the client.
COMPRESS_LEVEL = 6

# storing consensus tree file spit out by sum_trees.py
TEMP_CONSENSUS_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'temp_consensus')
CONSENSUS_PARALLIZATION = '4'

# for uploaded datasets
DATA_FILES_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'raw_data_files')
TREE_COLLECTION_FILENAME = 'tree_collection.tre'
TREE_COLLECTION_NAMES_FILENAME = 'tree_collection_names.txt'

# celery settings (same for dev and prod)
CELERY_BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'

# LOGIN_DISABLED = True