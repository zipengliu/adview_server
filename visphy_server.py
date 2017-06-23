from flask import Flask, jsonify
from flask_compress import Compress
from flask_cors import CORS
import sys
from pymongo import MongoClient


class Connection():
    def __init__(self, env):
        self.url = 'mongodb://visphybe.cs.ubc.ca:27017/visphy_dev' if env == 'prod' else 'mongodb://localhost:27017/visphy_dev'
        self.mongo_client = MongoClient(self.url)
        self.mongo_db = self.mongo_client.visphy_dev
        self.input_group = self.mongo_db.inputGroup
        self.entity = self.mongo_db.entity
        self.tree = self.mongo_db.tree
        self.branch = self.mongo_db.branch

    def test_connection(self):
        n = self.input_group.find({}).count()
        print 'MongoDB connected.'
        print 'There are', n, 'datasets in the database.'


env = 'dev'
if len(sys.argv) > 1:
    env = sys.argv[1]
connection = Connection(env)
connection.test_connection()
app = Flask('Visphy Server ' + env)
Compress(app)
CORS(app)


@app.route('/')
def hello_world():
    return 'Hello World!'


@app.route('/datasets')
def get_datasets():
    cursor = connection.input_group.find({})
    data = [{'inputGroupId': d['inputGroupId'], 'title': d['title'], 'description': d.get('description', ''),
             'numTrees': len(d['trees'])}
            for d in cursor]
    return jsonify(data)


@app.route('/dataset/<int:input_group_id>')
def get_dataset(input_group_id):
    print 'Getting dataset', input_group_id
    data = connection.input_group.find_one({'inputGroupId': input_group_id}, projection={'trees': False})
    entity_cursor = connection.entity.find({'inputGroupId': input_group_id}, projection={'eid': True, 'name': True})
    tree_cursor = connection.tree.find({'inputGroupId': input_group_id},
                                       projection={'name': True, 'tid': True, 'entities': True, 'rfDistance': True, 'rootBranch': True})
    branch_cursor = connection.branch.find({'inputGroupId': input_group_id},
                                           projection={'inputGroupId': False, 'cb': False, 'cb2': False, 'parent': False, 'isLeaf': False})
    ref_branch_cursor = connection.branch.find({'inputGroupId': input_group_id, 'tid': data['defaultReferenceTree']},
                                               projection={'inputGroupId': False, 'tid': False, 'parent': False, 'isLeaf': False})

    del data['_id']
    trees = {}
    for d in tree_cursor:
        trees[d['tid']] = d
        trees[d['tid']]['branches'] = {}
        del d['_id']
    for d in branch_cursor:
        trees[d['tid']]['branches'][d['bid']] = d
        del d['_id']
        del d['bid']
        del d['tid']
    entities = {}
    for d in entity_cursor:
        entities[d['eid']] = d
        del d['_id']

    ref_tree = trees.pop(data['defaultReferenceTree'])
    for d in ref_branch_cursor:
        ref_tree['branches'][d['bid']] = d
        del d['_id']
        del d['bid']

    data.update({
        'trees': trees,
        'entities': entities,
        'referenceTree': ref_tree
    })

    return jsonify(data)

if __name__ == '__main__':
    app.run(port=33333)
