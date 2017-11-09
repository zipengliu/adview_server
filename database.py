import config
from pymongo import MongoClient
from pymongo.errors import DocumentTooLarge
from gridfs import GridFS
import json


class Connection():
    def __init__(self, env, connect=True):
        self.url = config.MONGODB_URL_PROD if env == 'production' else config.MONGODB_URL_DEV
        print 'Connecting to MongoDB: ', self.url
        self.mongo_client = MongoClient(self.url, connect=connect)
        self.db_name = 'visphy_prod'if env == 'production' else 'visphy_dev'
        self.mongo_db = self.mongo_client[self.db_name]
        self.input_group = self.mongo_db.inputGroup
        self.entity = self.mongo_db.entity
        self.tree = self.mongo_db.tree
        self.branch = self.mongo_db.branch
        self.user = self.mongo_db.user
        self.gfs = GridFS(self.mongo_db)

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

    # Use GridFS to store document with size > 16MB
    def insert_gridfs(self, data):
        new_file = self.gfs.new_file()
        new_file.write(json.dumps(data))
        new_file.close()
        return new_file

    def read_grid_file(self, id):
        f = self.gfs.get(id)
        return json.load(f)

    def insert_tree(self, data):
        try:
            self.tree.insert_one(data)
        except DocumentTooLarge:
            tree_file = self.insert_gridfs(data)

            # Store some metadata with normal document
            self.tree.insert_one({'inputGroupId': data['inputGroupId'], 'tid': data['tid'], 'gfsFileId': tree_file._id})

    def insert_branches(self, branches):
        try:
            self.branch.insert_many(branches)
        except DocumentTooLarge:
            for b in branches:
                branch_file = self.insert_gridfs(b)
                self.branch.insert_one({'inputGroupId': b['inputGroupId'], 'tid': b['tid'], 'bid': b['bid'],
                                        'gfsFileId': branch_file._id})
