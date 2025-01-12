"""=================================================================================================
Diagon

Flask application for dotplots.
TODO sequence reverse currently acts like a toggle, if you make one reverse plot, and
then another, the sequence is forward in the second


Michael Gribskov     16 July 2020
================================================================================================="""
import sys
import os
import copy

from diagonal import Diagonal
from sequence.fasta import Fasta

from flask import Flask, render_template, request, redirect, flash, url_for, session
from flask_session import Session
from cachelib.file import FileSystemCache
import logging

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE

# from bokeh.util.string import encode_utf8

cli = sys.modules['flask.cli']
cli.show_server_banner = lambda *x: None

app = Flask(__name__)
app.secret_key = 'secret'
SESSION_TYPE = 'cachelib'
SESSION_SERIALIZATION_FORMAT = 'json'
SESSION_CACHELIB = FileSystemCache(threshold=500, cache_dir=f"{os.path.dirname(__file__)}/sessions")
USE_PERMANENT_SESSION = False

app.config.from_object(__name__)
Session(app)


# uncomment below to turn off server log
# log = logging.getLogger('werkzeug')
# log.setLevel(logging.ERROR)

# state is used to pass information to templates
def tf(string):
    if string == 'True':
        return True
    else:
        return False


statedefault = {'seq':     [{'fasta': None, 'status': 'next'},
                            {'fasta': None, 'status': 'later'}],
                'error':   '',
                'DNA':     {'identity': 'table/NUCidentity.matrix',
                            'NUC4.4':   'table/NUC4.4.matrix'},
                'protein': {'identity': 'table/PROidentity.matrix',
                            'Blosum62': 'table/BLOSUM62.matrix'},
                'params':  {'advanced':   False,
                            'mindotsize': 2,
                            'maxdotsize': 4,
                            'mode':       'dot',
                            'random':     'True',
                            'cmp':        'identity',
                            'width':      1,
                            'color':      1,
                            'window':     20,
                            'threshold':  12,
                            'seqtype':    'DNA',
                            'plottype':   'forward',
                            'cbase':      'Viridis'}
                }

state = copy.deepcopy(statedefault)


def getParams(res, state):
    for v in request.form:
        state['params'][v] = request.form[v]

    return


def isACGT(fasta, threshold=0.8):
    """-----------------------------------------------------------------------------------------
    Return True if at least threshold fraction of characters in the sequence are ACGT

    :param fasta: string    the sequence and nothing else
    :param threshold: float
    :return: Boolean
    -----------------------------------------------------------------------------------------"""
    total = len(fasta['seq'])
    if not total:
        return False

    # get the composition, assum uppercase
    count = {}
    for ch in fasta['seq']:
        if ch in count:
            count[ch] += 1
        else:
            count[ch] = 1

    # figure out fraction that are ACGT
    acgt = 0
    for base in 'ACGT':
        try:
            acgt += count[base]
        except KeyError:
            # ignore missing bases (but they count as non-ACGT)
            continue

    if acgt / total >= threshold:
        return True

    return False


# --------------------------------------------------------------------------------------------------
# Flask routes
# --------------------------------------------------------------------------------------------------

@app.route('/')
def index():
    """---------------------------------------------------------------------------------------------
    Generate a random session key and pass to dashboard.html

    :return: render dashboard.html
    ---------------------------------------------------------------------------------------------"""
    session_key = str(os.urandom(12))
    session[session_key] = {'state': copy.deepcopy(statedefault)}
    return render_template('dashboard.html', user=session_key)


@app.route('/dashboard', methods=['POST', 'GET'])
def dashboard():
    """---------------------------------------------------------------------------------------------
        present the dashboard: covers all stages of sequnce entry and parameter selection

        :return: render dashboard.html
    ---------------------------------------------------------------------------------------------"""
    return render_template('dashboard.html')


@app.route('/advanced', methods=['POST', 'GET'])
def advanced():
    # ----------------------------------------------------------------------------------------------
    # advanced toggles showing of the advanced parameters and reloads the dashboard
    # ----------------------------------------------------------------------------------------------

    sid = request.form.get("session_key")
    state = session[sid]['state']
    getParams(request, state)
    if state['params']['advanced']:
        state['params']['advanced'] = False
    else:
        state['params']['advanced'] = True

    return render_template('dashboard.html', state=state)


@app.route('/getSequence', methods=['POST', 'GET'])
def getSequence():
    # ----------------------------------------------------------------------------------------------
    # load sequences, store in session
    # ----------------------------------------------------------------------------------------------

    sequence = None
    if request.method == 'POST':
        sid = request.form.get("session_key")
        state = session[sid]['state']
        if 'file1' in request.files:
            f = request.files['file1']
            sequence = 0
            state['seq'][1]['status'] = 'next'

        elif 'file2' in request.files:
            f = request.files['file2']
            sequence = 1

        fasta = Fasta(fh=f)
        success = fasta.read()
        if not success:
            app.logger.warning('Cannot read file as fasta sequence')
            state['error'] = 'not a valid sequence file'
            return redirect(url_for('getSequence'))

        state['error'] = ''

        # store the fasta information in the session. unfortunately the object itself cannot be stored
        # sequence is 0 or 1 for sequnece 1 and 2, respectively
        print(fasta.format())
        seq = state['seq'][sequence]
        seq['fasta'] = {'id': fasta.id, 'doc': fasta.doc, 'seq': fasta.seq, 'format': fasta.format(), 'dir': 'forward'}
        seq['status'] = 'loaded'

        # else:
        #     # flash('Error in sequence file')
        #     # return render_template('dashboard.html', state=state)
        #     print(f'f:{f}\tlen:{f.content_length}\ttype:{f.content_type}\tparams'
        #                                  f':{f.mimetype_params}')
        #     state['error'] = 'not a valid sequence file'
        #     return redirect(url_for('getSequence'))

        # if both sequences have been selected, check whether the sequences are DNA or protein
        state['params']['seqtype'] = 'protein'
        state['params']['cmp'] = 'Blosum62'
        if state['seq'][0]['status'] == 'loaded' and state['seq'][1]['status'] == 'loaded':
            if isACGT(state['seq'][0]['fasta']) and isACGT(state['seq'][1]['fasta']):
                state['params']['seqtype'] = 'DNA'
                state['params']['cmp'] = 'identity'

        session.modified = True

    return render_template('dashboard.html')


@app.route('/self', methods=['POST', 'GET'])
def self():
    """---------------------------------------------------------------------------------------------
    when self dotplot is selected for the second sequence, this copies the first sequence into the
    second and returns to the main dashboard
    for a self plot, sequence 1 should already be loaded, we just have to copy it to sequence 2.
    :return:
    ---------------------------------------------------------------------------------------------"""
    sid = request.form.get("session_key")
    state = session[sid]['state']

    seq1 = state['seq'][0]
    seq2 = state['seq'][1]
    seq2['fasta'] = seq1['fasta'].copy()
    seq2['status'] = 'loaded'

    state['seqtype'] = 'protein'
    if state['seq'][0]['status'] == 'loaded' and state['seq'][1]['status'] == 'loaded':
        if isACGT(state['seq'][0]['fasta']) and isACGT(state['seq'][1]['fasta']):
            state['seqtype'] = 'DNA'

    session.modified = True

    return render_template('dashboard.html', state=state)


# --------------------------------------------------------------------------------------------------
# create the plot
# --------------------------------------------------------------------------------------------------
@app.route('/dotplot', methods=['POST', 'GET'])
def dotplot():
    sid = request.form.get("session_key")
    state = session[sid]['state']

    getParams(request, state)
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

    # int() cannot cast a not integer string, so must use int(float('1.5'))
    match.setupCalculation(fasta1, fasta2,
                           window=int(p['window']), threshold=int(float(p['threshold'])))
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
                               window=int(p['window']), threshold=int(float(p['threshold'])))
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

        # match.writeFrame('scoredist', key='score')
        # match.writeFrame('rundist', key='len')

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
    app.run(debug=True)

    exit(0)
