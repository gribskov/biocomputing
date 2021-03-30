"""=================================================================================================
Result example from Flask tutorial (modified)

Michael Gribskov     30 March 2021
================================================================================================="""
from flask import Flask, render_template

app = Flask(__name__)


@app.route('/result')
def result():
    dict = {'physics': 50, 'chemisty': 60, 'math': 70, 'python':97}
    return render_template('result.html', result=dict)


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run(debug=True)

    exit(0)
