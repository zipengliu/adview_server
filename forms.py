from wtforms import Form, StringField, PasswordField, BooleanField, FileField
from wtforms.validators import DataRequired, AnyOf


class SignupForm(Form):
    username = StringField('username', validators=[DataRequired()])
    password = PasswordField('password', validators=[DataRequired()])


class LoginForm(Form):
    username = StringField('username', validators=[DataRequired()])
    password = PasswordField('password', validators=[DataRequired()])


class DatasetForm(Form):
    title = StringField('title', validators=[DataRequired()])
    description = StringField('description', validators=[DataRequired()])
    support_values = StringField('supportValues', validators=[DataRequired()])
    is_reference_rooted = BooleanField('isReferenceRooted', false_values='N', validators=[DataRequired(), AnyOf(['Y', 'N'])])
    is_tc_rooted = BooleanField('isTCRooted', false_values='N', validators=[DataRequired(), AnyOf(['Y', 'N'])])
    reference_tree_file = FileField('reference', validators=[DataRequired()])
    tree_collection_file = FileField('treeCollection', validators=[DataRequired()])
    tree_collection_name_file = FileField('treeCollectionNames', validators=[])
    is_public = BooleanField('isPublic', false_values='N', validators=[DataRequired(), AnyOf(['Y', 'N'])])
