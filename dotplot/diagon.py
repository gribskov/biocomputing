"""=================================================================================================
Diagon

Flask application for dotplots.
TODO sequence reverse currently acts like a toggle, if you make one reverse plot, and
then another, the sequence is forward in the second


Michael Gribskov     16 July 2020
================================================================================================="""
import sys
import copy

from diagonal import Diagonal
from sequence.fasta import Fasta

from flask import Flask, render_template, request, redirect
import logging

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE

# from bokeh.util.string import encode_utf8

cli = sys.modules['flask.cli']
cli.show_server_banner = lambda *x: None
app = Flask(__name__)


# uncomment below to turn off server log
# log = logging.getLogger('werkzeug')
# log.setLevel(logging.ERROR)

# state is used to pass information to templates
def tf(string):
    if string == 'True':
        return True
    else:
        return False


statedefault = {'seq': [{'fasta': None, 'status': 'next'},
                        {'fasta': None, 'status': 'later'}],
                'DNA': {'identity': 'table/NUCidentity.matrix',
                           'NUC4.4': 'table/NUC4.4.matrix'},
                'protein': {'identity': 'table/PROidentity.matrix',
                           'Blosum62': 'table/BLOSUM62.matrix'},
                'params': {'advanced': False,
                           'mindotsize': 2,
                           'maxdotsize': 8,
                           'mode': 'dot',
                           'random': 'True',
                           'cmp': 'identity',
                           'width': 1,
                           'color': 1,
                           'window': 20,
                           'threshold': 12,
                           'seqtype': 'DNA',
                           'plottype': 'forward',
                           'cbase': 'Viridis'}
                }

state = copy.deepcopy(statedefault)


def getParams(res):
    for v in request.form:
        state['params'][v] = request.form[v]

    return


# --------------------------------------------------------------------------------------------------
# Flask routes
# --------------------------------------------------------------------------------------------------

@app.route('/')
def index():
    return redirect('/dashboard', code=302)


# --------------------------------------------------------------------------------------------------

@app.route('/dashboard', methods=['POST', 'GET'])
def dashboard():
    state = copy.deepcopy(statedefault)
    return render_template('dashboard.html', state=state)


# --------------------------------------------------------------------------------------------------
# advanced toggles showing of the advanced parameters and reloads the dashboard
# --------------------------------------------------------------------------------------------------
@app.route('/advanced', methods=['POST', 'GET'])
def advanced():
    getParams(request)
    if state['params']['advanced']:
        state['params']['advanced'] = False
    else:
        state['params']['advanced'] = True

    return render_template('dashboard.html', state=state)


# --------------------------------------------------------------------------------------------------
# main dashboard for loading sequences, scoring table and setting parameters
# --------------------------------------------------------------------------------------------------
@app.route('/getSequence', methods=['POST', 'GET'])
def getSequence():
    sequence = None
    if request.method == 'POST':
        if 'file1' in request.files:
            f = request.files['file1']
            sequence = 0
            state['seq'][1]['status'] = 'next'

        elif 'file2' in request.files:
            f = request.files['file2']
            sequence = 1

        fasta = Fasta(fh=f)
        fasta.read()
        print(fasta.format())
        seq = state['seq'][sequence]
        seq['fasta'] = fasta
        seq['status'] = 'loaded'

        # if both sequences have been selected, check whether the sequences are DNA or protein
        state['params']['seqtype'] = 'protein'
        if state['seq'][0]['status'] is 'loaded' and state['seq'][1]['status'] is 'loaded':
            if state['seq'][0]['fasta'].isACGT() and state['seq'][1]['fasta'].isACGT():
                state['params']['seqtype'] = 'DNA'

    return render_template('dashboard.html', state=state)


# --------------------------------------------------------------------------------------------------
# when self dotplot is selected for the second sequence, this copies the first sequence into the
# second and returns to the main dashboard
# --------------------------------------------------------------------------------------------------
@app.route('/self', methods=['POST', 'GET'])
def self():
    """
    for a self plot, sequence 1 should already be loaded, we just have to copy it to sequence 2.
    :return:
    """
    seq1 = state['seq'][0]
    seq2 = state['seq'][1]
    seq2['fasta'] = seq1['fasta'].copy()
    seq2['status'] = 'loaded'

    state['seqtype'] = 'protein'
    if state['seq'][0]['status'] is 'loaded' and state['seq'][1]['status'] is 'loaded':
        if state['seq'][0]['fasta'].isACGT() and state['seq'][1]['fasta'].isACGT():
            state['seqtype'] = 'DNA'

    return render_template('dashboard.html', state=state)


# --------------------------------------------------------------------------------------------------
# create the plot
# --------------------------------------------------------------------------------------------------
@app.route('/dotplot', methods=['POST', 'GET'])
def dotplot():
    getParams(request)
    p = state['params']

    mode = p['mode']
    plottype = p['plottype']
    seqtype = p['seqtype']

    fasta1 = state['seq'][0]['fasta']
    fasta2 = state['seq'][1]['fasta']

    match = Diagonal()
    color = int(p['color'])
    width = int(p['width'])
    match.mindotsize = int(p['mindotsize'])
    match.maxdotsize = int(p['maxdotsize'])

    cmpname = p['cmp']
    cmptable = state[seqtype][cmpname]
    match.readNCBI(cmptable)

    dataframes = [{'data': 'dots', 'fn': match.windowThreshold, 'var': ['x', 'y', 'score']},
                  {'data': 'scoredist', 'fn': match.histogramScore, 'var': ['score', 'count']},
                  {'data': 'rundist', 'fn': match.histogramRun, 'var': ['len', 'count']},
                  {'data': 'randomscore', 'fn': None, 'var': ['score', 'count']},
                  {'data': 'randomrun', 'fn': None, 'var': ['len', 'count']}
                  ]
    match.setupFrame(dataframes)

    if plottype == "reverse":
        match.seqreverse = True
    match.setupCalculation(fasta1, fasta2,
                           window=int(p['window']), threshold=int(p['threshold']))
    match.setupBokeh(cbase=p['cbase'], clevels=256, creverse='True')
    match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
    if mode == 'line':
        match.addSegment('dots')
    match.bdot('dots', 'main', width=width, color=color,
               mode=p['mode'])

    if plottype == "forward_backward":
        match.seqreverse = True
        match.resetFrame('dots')
        match.setupCalculation(fasta1, fasta2, resetstat=False,
                               window=int(p['window']), threshold=int(p['threshold']))
        # match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
        match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
        if mode == 'line':
            match.addSegment('dots')
        match.bdot('dots', 'main', set_colormap=False,
                   width=width, color=color, mode=p['mode'])

    # score and run distributions
    match.sortFrame('scoredist', 'score')
    match.addCumulative('scoredist', 'count', 'cumulative')
    match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
    match.bscoreCumulative('scoredist', 'scoredist')
    match.brunDist('rundist', 'rundist', color='#0000ff')

    # random score distribution
    if p['random'] == 'True':
        match.single = True
        match.random(n=match.nscore)
        match.histogramScore('randomscore', 1)
        match.sortFrame('randomscore', 'score')
        match.bscoreDist('scoredist', 'randomscore', color='#ff0000')

        match.histogramRun('randomrun', 1)
        match.brunDist('rundist', 'randomrun', color='#ff0000')
        match.single = False

        match.writeFrame('scoredist', key='score')
        match.writeFrame('rundist', key='len')

    script, div = components(match.grid)
    # grab the static resources
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    return render_template(
        'dotplot.html',
        plot_script=script,
        plot_div=div,
        js_resources=js_resources,
        css_resources=css_resources
    )


# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # print('Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)')
    app.run(debug=True)
    # app.run()

    exit(0)
