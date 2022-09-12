"""=================================================================================================
Result example from Flask tutorial (modified)

Michael Gribskov     30 March 2021
================================================================================================="""
from flask import Flask, render_template, request
app = Flask(__name__)

@app.route('/')
def enter_scores():
   return render_template('enter_scores.html')

@app.route('/result',methods = ['POST', 'GET'])
def result():
   if request.method == 'POST':
      result = request.form
      return render_template("result.html",result = result)


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run(debug=True)

    exit(0)
