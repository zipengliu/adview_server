import sys, json
from pymongo import MongoClient

mongo_client = MongoClient('localhost', 27017)
mongo_db = mongo_client.visphy_dev
entity_col = mongo_db.entity

input_group_id = int(sys.argv[1])


mapping_filename = '../data/conifer_mapping.json'
m = json.load(open(mapping_filename))

cnt = 0
for k in m:
    r = entity_col.update_one({'inputGroupId': input_group_id, 'name': k}, {'$set': {'name': m[k]}})
    if r.modified_count: 
        cnt += 1
        print k, '->', m[k]


print cnt, 'entity names changed.'
print 'Done'
