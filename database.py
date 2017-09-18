import config
from pymongo import MongoClient


class Connection():
    def __init__(self, env):
        self.url = config.MONGODB_URL_PROD if env == 'production' else config.MONGODB_URL_DEV
        self.mongo_client = MongoClient(self.url)
        self.mongo_db = self.mongo_client.visphy_dev
        self.input_group = self.mongo_db.inputGroup
        self.entity = self.mongo_db.entity
        self.tree = self.mongo_db.tree
        self.branch = self.mongo_db.branch
        self.user = self.mongo_db.user

    def test_connection(self):
        n = self.input_group.find({}).count()
        print 'MongoDB connected.'
        print 'There are', n, 'datasets in the database.'

    # This function is prone to race condition
    def get_next_input_group_id(self):
        all_input_groups = self.input_group.find(projection={'inputGroupId': True})
        max_id = 0
        for d in all_input_groups:
            if d['inputGroupId'] > max_id:
                max_id = d['inputGroupId']
        return max_id + 1
