"""=================================================================================================
Example Flask application

Michael Gribskov     30 March 2021
================================================================================================="""
from flask import Flask, request, redirect, url_for
app = Flask(__name__)

@app.route('/<worldtype>')
def welcome(worldtype):
    return '<H1>Welcome to the world of {}<H1>'.format(worldtype)

@app.route('/login',methods = ['POST', 'GET'])
def login():
   if request.method == 'POST':
       user = request.form['user']
       return redirect(url_for('welcome', worldtype=user))
   else:
      user = request.args.get('user')
      url = '/{}'.format(user)
      return redirect(url)


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    app.debug = True
    app.run()

    exit(0)
