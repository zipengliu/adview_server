import config
import secret_key
from flask import Flask, jsonify, request, url_for
from flask_compress import Compress
from flask_cors import CORS
from flask_login import LoginManager, UserMixin, login_user, logout_user, login_required, current_user
from werkzeug.utils import secure_filename
from dendropy import Tree, TreeList
import os
import subprocess
from functools import wraps
from datetime import datetime
from passlib.context import CryptContext
from bson.objectid import ObjectId
import uuid

from database import Connection
from forms import SignupForm, LoginForm, DatasetForm
from preprocess_worker import preprocess_dataset


# Initialize the server application
app = Flask('Visphy Server')
app.config.from_object('config')
app.config.from_envvar('ENVIRONMENT', silent=True)
app.config.from_object('secret_key')
pwd_crpt = CryptContext(["sha256_crypt", "ldap_salted_md5"])

login_manager = LoginManager()
login_manager.init_app(app)


class User(UserMixin):

    def __init__(self, _id, username, password):
        self._id = _id
        self.username = username
        self.password = password

    @staticmethod
    def get_by_username(username):
        user = connection.user.find_one({'username': username})
        return User.get_from_document(user)

    @staticmethod
    def get_by_id(id):
        user = connection.user.find_one({'_id': ObjectId(id)})
        return User.get_from_document(user)

    @staticmethod
    def get_from_document(d):
        if d is None:
            return d
        return User(str(d['_id']), d['username'], d['password'])

    def save(self):
        connection.user.insert_one({'username': self.username, 'password': self.password})

    def get_id(self):
        return unicode(self._id)


@login_manager.user_loader
def load_user(user_id):
    return User.get_by_id(user_id)


@login_manager.unauthorized_handler
def unauthorized():
    return jsonify({'status': 'Unauthorized'}), 401


@app.route('/signup', methods=['POST'])
def register():
    form = SignupForm(request.form)
    if form.validate():
        # Check existent username
        if connection.user.find_one({'username': form.username.data}) is not None:
            return 'User already exist.', 400

        user = User('', form.username.data, pwd_crpt.hash(form.password.data))
        user.save()
        return "User successfully registered."
    else:
        return 'Invalid form', 400


@app.route('/login', methods=['POST'])
def login():
    form = LoginForm(request.form)
    if form.validate():
        user = User.get_by_username(form.username.data)
        if user:
            if pwd_crpt.verify(form.password.data, user.password):
                # successfully logged in
                login_user(user)
                return 'User logged in'
            else:
                # password not correct
                return 'Username or password is not correct', 400
        else:
            return 'Cannot find this user', 400
    else:
        return 'Invalid form', 400


# Decorator for authorizing a dataset
# Needs to be called after login_required()
def authorize_dataset(func):
    @wraps(func)
    def authorize_and_call(*args, **kwargs):
        if not app.config['LOGIN_DISABLED']:
            input_group_id = request.args.get('inputGroupId', None) or request.form.get('inputGroupId', None)
            if not input_group_id:
                return 'No dataset id provided', 400

            d = connection.input_group.find_one({'inputGroupId': int(input_group_id)},
                                                projection={'ownerUserId': True, 'isPublic': True})
            if d is None:
                return 'There is no such dataset!', 400
            if not d.get('isPublic', False) and str(d.get('ownerUserId', '')) != current_user._id:
                return 'Unauthorized access!', 403
        return func(*args, **kwargs)
    return authorize_and_call


@app.route('/logout')
@login_required
def logout():
    logout_user()
    return 'User logged out'


@app.route('/datasets')
@login_required
def get_datasets():
    fields = ['inputGroupId', 'title', 'description', 'numTrees', 'numTaxa', 'referenceTreeFileName', 'timeCreated', 'isPublic']

    def get_datasets_by_filter(f):
        cursor = connection.input_group.find(f)
        return [{f:(d.get(f, 'N/A') if not f.startswith('time') else d['_id'].generation_time) for f in fields}
                for d in cursor]
    data = get_datasets_by_filter({'ownerUserId': ObjectId(current_user._id)}) + get_datasets_by_filter({'isPublic': True})
    return jsonify(data)


@app.route('/dataset', methods=['GET'])
@authorize_dataset
@login_required
def get_dataset():
    input_group_id = request.args.get('inputGroupId', None)
    if input_group_id is None:
        return 'No dataset id is provided', 400
    input_group_id = int(input_group_id)

    data = connection.input_group.find_one({'inputGroupId': input_group_id}, projection={'trees': False, '_id': False, 'ownerUserId': False})
    if data is None:
        return 'No such dataset', 400

    entity_cursor = connection.entity.find({'inputGroupId': input_group_id}, projection={'eid': True, 'name': True, '_id': False})
    tree_cursor = connection.tree.find({'inputGroupId': input_group_id},
                                       projection={'name': True, 'tid': True, 'entities': True, 'rfDistance': True, 'rootBranch': True, '_id': False, 'outgroupBranch': True})
    branch_cursor = connection.branch.find({'inputGroupId': input_group_id},
                                           projection={'inputGroupId': False, 'cb': False, 'cb2': False, 'parent': False, 'isLeaf': False, '_id': False})
    ref_branch_cursor = connection.branch.find({'inputGroupId': input_group_id, 'tid': data.get('defaultReferenceTree', 't0')},
                                               projection={'inputGroupId': False, 'tid': False, 'parent': False, 'isLeaf': False, '_id': False})

    trees = {}
    for d in tree_cursor:
        trees[d['tid']] = d
        trees[d['tid']]['branches'] = {}
    for d in branch_cursor:
        trees[d['tid']]['branches'][d['bid']] = d
        del d['bid']
        del d['tid']
    entities = {}
    for d in entity_cursor:
        entities[d['eid']] = d

    ref_tree = trees.pop(data.get('defaultReferenceTree', 't0'))
    for d in ref_branch_cursor:
        ref_tree['branches'][d['bid']] = d
        del d['bid']

    data.update({
        'trees': trees,
        'entities': entities,
        'referenceTree': ref_tree
    })

    return jsonify(data)


@app.route('/dataset', methods=['DELETE'])
@authorize_dataset
@login_required
def delete_dataset():
    input_group_id = request.args.get('inputGroupId', None)
    if input_group_id is None:
        return 'No dataset id is provided', 400
    input_group_id = int(input_group_id)

    connection.entity.delete_many({'inputGroupId': input_group_id})
    connection.tree.delete_many({'inputGroupId': input_group_id})
    connection.branch.delete_many({'inputGroupId': input_group_id})
    connection.input_group.delete_one({'inputGroupId': input_group_id})
    return 'success'


def read_tree_from_consensus_report(filepath):
    newick_string = None
    with open(filepath, 'r') as f:
        mark = False
        for line in f:
            if mark:
                newick_string = line.strip()
                break
            if line.startswith('BEGIN TREES'):
                mark = True
    newick_string = newick_string[newick_string.find('('):]
    return newick_string


def parse_tree(tree, entities):
    tree.ladderize()

    branches = {}
    no = 0
    for node in tree.levelorder_node_iter():
        assert node.is_leaf() or node.num_child_nodes() == 2
        node.bid = 'b' + str(no)
        no += 1

    present_entities = {}
    for node in tree.levelorder_node_iter():
        if node.is_leaf():
            # label = node.taxon.label.replace(' ', '_')
            label = node.taxon.label
            branches[node.bid] = {
                'entity': entities[label],
                # 'support': node.support,
                'length': node.edge_length or 0
            }
            present_entities[entities[label]] = True
        else:
            support_val = 0 if not node.label else float(node.label)
            branches[node.bid] = {
                'left': node.child_nodes()[0].bid,
                'right': node.child_nodes()[1].bid,
                'support': support_val,
                'length': node.edge_length or 0
            }

    data = {
        'name': 'Consensus tree',
        'branches': branches,
        'entities': present_entities.keys(),
        'rfDistance': {},       # Cost is too high?
        'rootBranch': 'b0',
        'tid': 'consensus-' + str(uuid.uuid4())
    }

    return data


# Retrieve the tree data (newick string) from the DB
def get_newick_from_DB(input_group_id, tids):
    tree_cursor = connection.tree.find({'inputGroupId': input_group_id}, projection={'tid': True, 'newickString': True})
    trees = {tid: True for tid in tids}
    remain = len(tids)
    tree_str = []
    for d in tree_cursor:
        if d['tid'] in trees:
            tree_str.append(d['newickString'])
            remain -= 1
            if remain == 0:
                break
    assert remain == 0, 'Failed to find newick string for some trees'
    return tree_str


@app.route('/consensus', methods=['POST'])
@authorize_dataset
@login_required
def get_consensus():
    input_group_id = int(request.form['inputGroupId'])
    tids = str(request.form['trees']).split(',')
    print 'Creating consensus tree for ', tids, 'in dataset', input_group_id

    tree_str = get_newick_from_DB(input_group_id, tids)
    ### This should be as easy as to call the TreeList.consensus method, but I don't know why the summerized tree
    ### generated by this function does not have edge length.
    # tree_list = TreeList([Tree.get(data=s, schema='newick') for s in tree_str])
    # tree = tree_list.consensus(min_freq=0.6)
    # print tree.length()

    # Write the strings to a file
    tree_filename = datetime.now().strftime('%m%d%H%M%S%f') + '.tre'
    tree_path = os.path.join(app.config['TEMP_CONSENSUS_PATH'], tree_filename)
    with open(tree_path, 'w') as f:
        f.write('\n'.join(tree_str))

    # Invoke sumtrees.py
    consensus_path = os.path.join(app.config['TEMP_CONSENSUS_PATH'], tree_filename + '.consensus')
    try:
        subprocess.check_call(['sumtrees.py', '--output-tree-filepath=' + consensus_path, '--suppress-annotations',
                               '--multiprocessing=' + app.config['CONSENSUS_PARALLIZATION'], tree_path],
                              stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        return 'Error summerizing trees', 500

    # Read the entity mapping from DB
    entity_cursor = connection.entity.find({'inputGroupId': input_group_id}, projection={'eid': True, 'name': True})
    entities = {str(e['name']): str(e['eid']) for e in entity_cursor}
    #
    # # Read the consensus tree and parse
    newick_string = read_tree_from_consensus_report(consensus_path)
    tree = Tree.get(data=newick_string, schema='newick')

    tree.resolve_polytomies()
    parsed_tree = parse_tree(tree, entities)
    parsed_tree['newickString'] = tree.as_string('newick')
    return jsonify(parsed_tree)


@app.route('/tree_newick', methods=['POST'])
@authorize_dataset
@login_required
def get_tree_newick():
    input_group_id = int(request.form['inputGroupId'])
    tids = str(request.form['tids']).split(',')
    print 'Getting newick strings for ', tids, 'in dataset', input_group_id

    tree_str = get_newick_from_DB(input_group_id, tids)
    return ''.join(tree_str)


@app.route('/dataset', methods=['POST'])
@login_required
def create_new_dataset():
    print 'Creating new dataset...'
    print request.form
    # form = DatasetForm(request.POST)
    # if not form.validate():
    #     return 'Invalid form', 400
    try:
        title = request.form['title']
        description = request.form['description']
        support_value_range = request.form['supportValues']
        is_reference_rooted = request.form['isReferenceRooted'] == 'Y'
        is_tc_rooted = request.form['isTCRooted'] == 'Y'
        reference_tree_file = request.files['reference']
        reference_tree_filename = secure_filename(reference_tree_file.filename)
        tree_collection_file = request.files['treeCollection']
        tree_collection_name_file = request.files.get('treeCollectionNames', None)      # optional
        is_public = request.form.get('isPublic', 'N') == 'Y'
    except Exception as e:
        print 'Dataset is not acceptable', e
        return 'Dataset is not acceptable', 400

    # Process the dataset
    # 1. Add an inputGroup in DB
    input_group_id = connection.get_next_input_group_id()
    input_group_data = {
        'inputGroupId': input_group_id,
        'title': title,
        'description': description,
        'isReferenceRooted': is_reference_rooted,
        'isTCRooted': is_tc_rooted,
        'supportRange': None if support_value_range == 'NA' else ([0, 1] if support_value_range == '1' else [0, 100]),
        # 'withParalogs': False,
        # 'hasMissingTaxa': False,        # Update this field after parsing
        'trees': [],
        'referenceTreeFileName': reference_tree_filename,
        'isPublic': is_public,
        'ownerUserId': ObjectId(current_user._id) if not app.config['LOGIN_DISABLED'] else None
    }

    # 2. store files in filesystem
    dirname = os.path.join(app.config['DATA_FILES_PATH'], str(input_group_id))
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    reference_tree_file.save(os.path.join(dirname, reference_tree_filename))
    tree_collection_file.save(os.path.join(dirname, app.config['TREE_COLLECTION_FILENAME']))
    if tree_collection_name_file:
        tree_collection_name_file.save(os.path.join(dirname, app.config['TREE_COLLECTION_NAMES_FILENAME']))
        tree_collection_name_file.close()

    # 3. Parse trees
    reference_tree_file.seek(0)
    reference_tree = Tree.get(data=reference_tree_file.read(), schema='newick')
    reference_tree.resolve_polytomies()
    tn = reference_tree.taxon_namespace
    print tn
    tree_collection_file.seek(0)
    tree_collection = TreeList.get(data=tree_collection_file.read(), schema='newick', taxon_namespace=tn)
    input_group_data['numTrees'] = len(tree_collection)
    input_group_data['numTaxa'] = len(tn)
    connection.input_group.insert_one(input_group_data)

    # 4. Close the stored files
    reference_tree_file.close()
    tree_collection_file.close()

    # 5. Construct the tree data structure for frontend display
    branches = {}
    no = 0
    for node in reference_tree.levelorder_node_iter():
        bid = 'b' + str(no)
        no += 1
        node.label = bid       # a quick hack

    for node in reference_tree.levelorder_node_iter():
        if node.is_internal():
            branches[node.label] = {
                'left': node.child_nodes()[0].label,
                'right': node.child_nodes()[1].label,
                'length': node.edge_length or 1
            }
        else:
            branches[node.label] = {
                'entity': 'e' + str(tn.accession_index(node.taxon)),
                'length': node.edge_length or 1
            }

    data = {
        'inputGroupId': input_group_id,
        'referenceTree': {
            'name': reference_tree_filename,
            'tid': 't0',
            'branches': branches
        },
        'entities': {'e' + str(tn.accession_index(t)): t.label for t in tn}
    }
    return jsonify(data)


# Hand the heavy lifting of preprocessing to the celery worker.
# Return a task ID for client to check status
@app.route('/outgroup', methods=['POST'])
@authorize_dataset
@login_required
def upload_outgroup():
    print 'Received outgroup data'
    input_group_id = int(request.form['inputGroupId'])
    is_updating = request.form.get('isUpdating', False)
    task = preprocess_dataset.delay(input_group_id, request.form['outgroup'], is_updating)
    print task
    return jsonify({'taskId': task.id}), 202, {'Location': url_for('check_upload_status', task_id=task.id)}


@app.route('/upload_status/<task_id>')
def check_upload_status(task_id):
    task = preprocess_dataset.AsyncResult(task_id)
    print task.info
    if task.state == 'PENDING':
        response = {
            'state': task.state,
        }
    elif task.state == 'SUCCESS':
        response = {
            'state': task.state,
            'url': '/dataset/' + str(task.info['input_group_id'])      # The view url for dataset
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
        }
        response.update(task.info)
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'error': str(task.info),
        }
    return jsonify(response)


connection = Connection(app.config['ENVIRONMENT'])
connection.test_connection()
Compress(app)
CORS(app, supports_credentials=True)

app.run(port=33333)
