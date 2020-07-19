"""=================================================================================================
Diagon

Flask application for dotplots.

Michael Gribskov     16 July 2020
================================================================================================="""
import sys

from diagonal import Diagonal
from sequence.fasta import Fasta

from flask import Flask, render_template, request
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE

# from bokeh.util.string import encode_utf8

cli = sys.modules['flask.cli']
cli.show_server_banner = lambda *x: None
app = Flask(__name__)

state = {'s1': {'fasta': None, 'isloaded': False},
         's2': {'fasta': None, 'isloaded': False}}


# ---------------------------------------------------------------------------------------------------
# Flask routes
# ---------------------------------------------------------------------------------------------------

@app.route('/')
def index():
    return render_template('dashboard.html', state=state)
    # return '<h1><a href=/bokeh>bokeh test</a></h1><a href=/dashboard>dashboard</a>'


@app.route('/dashboard')
def dashboard():
    return render_template('dashboard.html', state=state)


@app.route('/getSequence', methods=['POST', 'GET'])
def getSequence():
    sequence = None
    if request.method == 'POST':
        if 'file1' in request.files:
            f = request.files['file1']
            sequence = 's1'


        elif 'file2' in request.files:
            f = request.files['file2']
            sequence = 's2'

        fasta = Fasta(fh=f)
        fasta.read()
        print(fasta.format())
        state[sequence]['fasta'] = fasta
        state[sequence]['isloaded'] =  True

    return render_template('dashboard.html', state=state)


@app.route('/dotplot', methods=['POST', 'GET'])
def dotplot():
    match = Diagonal()

    print('testing')
    return 'testing'


@app.route('/bokeh')
def bokeh():
    # init a basic bar chart:
    # http://bokeh.pydata.org/en/latest/docs/user_guide/plotting.html#bars
    fig = figure(plot_width=600, plot_height=600)
    fig.vbar(
        x=[1, 2, 3, 4],
        width=0.5,
        bottom=0,
        top=[1.7, 2.2, 4.6, 3.9],
        color='navy'
    )

    # grab the static resources
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    # render template
    script, div = components(fig)
    html = render_template(
        'index.html',
        plot_script=script,
        plot_div=div,
        js_resources=js_resources,
        css_resources=css_resources,
    )
    return html


# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # state['seq1'] = None

    app.run(debug=True)
    # app.run()

    exit(0)
