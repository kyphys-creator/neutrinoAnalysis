#!/usr/bin/env python3
"""Restructure the auto-converted legacy .wl files into the QuantumSensor-style
numbered layout (NN_topic.wl, NNb_..._plot.wl), per threshold.

Input : <thr>/src/*.wl as written by convert_nb_to_wl.wls (with section markers)
Output: <thr>/src/00_main.wl ... 16_testing.wl  (legacy-named files removed)

Run from the Mathematica/ folder:  python3 tools/restructure_src.py
"""
import re
import pathlib

THRESHOLDS = ['1eV', '5eV']

HDR = ('(* {eq} *)\n'
       '(*  {name}  ({thr})\n*)\n'
       '(*  {prov}\n*)\n'
       '(* {eq} *)\n\n')
EQ = '=' * 74


def header(name, thr, prov):
    return HDR.format(eq=EQ, name=name, thr=thr, prov=prov)


def body_of(path):
    """File content minus the auto-generated header block and frontend calls."""
    text = path.read_text()
    # drop the 6-line converter header
    m = re.match(r'(?:\(\*.*?\*\)\n){1,6}\n', text, re.S)
    if m and 'Auto-converted' in m.group(0):
        text = text[m.end():]
    text = re.sub(r'\n?NotebookClose\[\]\n?', '\n', text)
    return text.strip() + '\n'


def blocks_of(text):
    """Split on blank lines into statement blocks (marker comments included)."""
    return [b for b in re.split(r'\n\s*\n', text) if b.strip()]


DIRS_BLOCK = '''(* ---- directories: src (this file), ../input, ../output ---- *)
Direc = DirectoryName[$InputFileName];
inputDir = FileNameJoin[{ParentDirectory[Direc], "input"}];
outputDir = FileNameJoin[{ParentDirectory[Direc], "output"}];
Quiet[CreateDirectory[FileNameJoin[{outputDir, #}], CreateIntermediateDirectories -> True]] & /@
  {"CR", FileNameJoin[{"CR_revised", "originalUnit"}]};
'''


def route_calelse(block):
    """Route one 5_Calculation block to a target file key."""
    first = block.lstrip()
    rules = [
        (r'^(datEvsolid|datdPhisolid|dfluxsolid\s*=|dfluxsolid$|intfluxsolid\s*=|'
         r'datEvdash|datdPhidash|dfluxdash\s*=|intfluxdash\s*=)', '06'),
        (r'^Plot\[\{(dfluxsolid|intfluxsolid)', '06b'),
        (r'^(ERbin2?\s*=|Ratebin7\s*=|Ratebin2\s*=|RatebinDif\s*=)', '07'),
        (r'^(Bin7plot\s*=|Bin2plot\s*=|BinDifplot\s*=|Show\[Bin7plot, Bin2plot)', '07b'),
        (r'^(CRbin7NoG|CRlistbox|CRlist\s*=|CRplotlist|CRbin7Box)', '08'),
        (r'^(Bin7plotCRNoG|Plot\[CRplotlist|Bin7plotCRBox|Show\[Bin7plot, Bin7plotCRNoG)', '08b'),
    ]
    for pat, key in rules:
        if re.match(pat, first):
            return key
    raise SystemExit(f'unrouted 5_Calculation block:\n{block[:120]}')


def split_main(text):
    """Split the (fixed) 1_Main body into named regions using section markers
    and content anchors. Returns dict of region-name -> text."""
    lines = text.splitlines(keepends=True)

    def find(pat, start=0):
        for i in range(start, len(lines)):
            if re.search(pat, lines[i]):
                return i
        raise SystemExit(f'anchor not found: {pat}')

    i_calc   = find(r'===== Section: Calculations =====')
    i_stored = find(r'===== Subsection: Stored Data =====')
    i_plots  = find(r'===== Section: Plots =====')
    i_curly  = find(r'===== Subsection: Curly R Matrix =====')
    i_sample = find(r'^SampleIntByEdges\[f_\]')
    i_n68    = find(r'===== Section: n=68 case =====')
    i_n100   = find(r'===== Section: n=100 case =====')
    i_n160   = find(r'===== Section: n=160 case =====')
    i_n160b  = find(r'===== Section: n=160 case with bkg')

    setting = ''.join(lines[:i_calc])
    checks = re.findall(r'^\(2\*\([^\n]+\n', setting, re.M)

    def rng(a, b):
        return ''.join(lines[a:b]).strip() + '\n'

    # Stored-Data leftovers: everything in the subsection except the (already
    # commented) Stored_Data note.
    stored = rng(i_stored, i_plots)
    stored = re.sub(r'\(\* ===== Subsection: Stored Data ===== \*\)\n', '', stored)
    stored = re.sub(r'\(\* Stored_Data\.nb.*?\*\)\n', '', stored, flags=re.S)

    return {
        'checks':    checks,
        'main_plots': stored.strip() + '\n\n' + rng(i_plots, i_curly),
        'ratebins_export': rng(i_curly, i_sample),
        'crmat_helpers':   rng(i_sample, i_n68),
        'n68':  rng(i_n68, i_n100),
        'n100': rng(i_n100, i_n160),
        'n160': rng(i_n160, i_n160b),
        'n160b': rng(i_n160b, len(lines)),
    }


MAIN_TMPL = '''srcDir = DirectoryName[$InputFileName];
load[f_] := Get[FileNameJoin[{{srcDir, f}}]];

(* ---- definitions & data (fast) ---- *)
load["01_setup.wl"];              (* directories, GeV scale, units *)
load["02_functions.wl"];          (* plot styles, smearing, integration helpers *)
load["03_constants_nuflux.wl"];   (* physics constants for the neutrino flux *)
load["04_exprate_defs.wl"];       (* expected-rate response functions *)
load["05_input_data.wl"];         (* recoil bins + dN/dE spectra from ../input *)

(* threshold checks (from the legacy Main "Setting" section) *)
{checks}
(* ---- flux, rate bins, curly-R (moderate NIntegrate load) ---- *)
load["06_flux_defs.wl"];
load["07_ratebins.wl"];
load["08_curlyR_defs.wl"];

(* ---- pipeline outputs -> ../output ---- *)
load["10_ratebins_export.wl"];    (* Ratebin7/2 (+ originalUnit) CSVs *)
load["11_crmat_helpers.wl"];      (* SampleIntByEdges / makeEvEdgesByCount *)

(* ---- optional: plots ---- *)
(* load["06b_flux_plot.wl"]; *)
(* load["07b_ratebins_plot.wl"]; *)
(* load["08b_curlyR_plot.wl"]; *)
(* load["09_main_plots.wl"]; *)

(* ---- optional: CRmat construction + chi^2 minimization studies ---- *)
(* load["12_minimization_n68.wl"]; *)
(* load["13_minimization_n100.wl"]; *)
(* load["14_minimization_n160.wl"]; *)
(* load["15_minimization_n160_bkg2MeV.wl"]; *)

(* ---- optional: consistency tests ---- *)
(* load["16_testing.wl"]; *)

(* Stored_Data.nb is not part of the legacy folder; re-enable if added:
   load["Stored_Data.wl"]  *)
'''


def process(thr):
    src = pathlib.Path(thr) / 'src'
    P = lambda n: src / n

    # ---------- component files (with path fixes) ----------
    units = body_of(P('Units.wl')).replace(
        '*"\\!\\(\\*SuperscriptBox[\\(cm\\), \\(-3\\)]\\)"', '*cm^(-3)')
    setup = (header('01_setup', thr,
                    'Directories, energy scale and units. From Units.nb; '
                    'GeV scale from 1_Main.nb (Setting).')
             + DIRS_BLOCK + '\n(* ---- energy scale (from 1_Main) ---- *)\n'
             + 'GeV = 10^15;\n\n' + units)

    functions = header('02_functions', thr, 'From Functions.nb.') \
        + body_of(P('Functions.wl'))
    constants = header('03_constants_nuflux', thr, 'From Constants_NuFlux.nb.') \
        + body_of(P('Constants_NuFlux.wl'))
    exprate = header('04_exprate_defs', thr, 'From ExpRate.nb.') \
        + body_of(P('ExpRate.wl'))

    inp = body_of(P('2_3_Input_Data.wl'))
    inp = inp.replace(
        'datdNdEsolid = Import[FileNameJoin[Direc, "dNdEsolid.csv"], "Data"];',
        'inputDir = FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]], "input"}];\n'
        'datdNdEsolid = Import[FileNameJoin[{inputDir, "dNdEsolid.csv"}], "Data"];')
    inp = inp.replace(
        'datdNdEdash = Import[FileNameJoin[Direc, "dNdEdash.csv"], "Data"];',
        'datdNdEdash = Import[FileNameJoin[{inputDir, "dNdEdash.csv"}], "Data"];')
    input_data = header('05_input_data', thr,
                        'Recoil-energy bins and dN/dE spectra. From 2_3_Input_Data.nb; '
                        'CSVs read from ../input.') + inp

    # ---------- split 5_Calculation ----------
    calc_targets = {k: [] for k in ('06', '06b', '07', '07b', '08', '08b')}
    for block in blocks_of(body_of(P('5_Calculation.wl'))):
        calc_targets[route_calelse(block)].append(block)
    calc_files = {
        '06_flux_defs.wl': ('Reactor-antineutrino flux interpolations '
                            '(dfluxsolid/dash, intfluxsolid/dash). From 5_Calculation.nb.', '06'),
        '06b_flux_plot.wl': ('Flux plots. From 5_Calculation.nb.', '06b'),
        '07_ratebins.wl': ('Recoil bins and expected-rate bins Ratebin7/2/Dif. '
                           'From 5_Calculation.nb.', '07'),
        '07b_ratebins_plot.wl': ('Rate-bin plots. From 5_Calculation.nb.', '07b'),
        '08_curlyR_defs.wl': ('Curly-R response: CRbin7NoG, CRloop/CRloopBox lists, '
                              'CRbin7Box. From 5_Calculation.nb.', '08'),
        '08b_curlyR_plot.wl': ('Curly-R plots. From 5_Calculation.nb.', '08b'),
    }

    # ---------- split 1_Main ----------
    main_body = body_of(P('1_Main.wl'))
    # legacy get-chain triplets & bare display lines in Setting are dropped by
    # region selection; rewrite exports to ../output first
    main_body = re.sub(r'FileNameJoin\[Direc, "test/([^"]+)"\]',
                       r'FileNameJoin[{outputDir, "\1"}]', main_body)
    main_body = re.sub(
        r'\w+ = NotebookOpen\[FileNameJoin\[Direc, "Stored_Data\.nb"\]\]\n'
        r'SelectionMove\[\w+, All, Notebook\]\nSelectionEvaluate\[\w+\]',
        '(* Stored_Data.nb is not part of the legacy folder *)', main_body)
    regions = split_main(main_body)

    main_files = {
        '09_main_plots.wl': ('Overview plots (flux data, expected rate, curly-R '
                             'matrix checks). From 1_Main.nb: "Stored Data" + "Plots".',
                             regions['main_plots']),
        '10_ratebins_export.wl': ('Ratebin7/2 at PrecisionGoal->4 and CSV exports '
                                  'to ../output. From 1_Main.nb: "Curly R Matrix".',
                                  regions['ratebins_export']),
        '11_crmat_helpers.wl': ('SampleIntByEdges / makeEvEdgesByCount helpers and '
                                'grid probes. From 1_Main.nb.',
                                regions['crmat_helpers']),
        '12_minimization_n68.wl': ('CRmat construction + chi^2 minimization, n=68 '
                                   'study. From 1_Main.nb: "n=68 case".', regions['n68']),
        '13_minimization_n100.wl': ('Same, n=100 study. From 1_Main.nb.', regions['n100']),
        '14_minimization_n160.wl': ('Same, n=160 study. From 1_Main.nb.', regions['n160']),
        '15_minimization_n160_bkg2MeV.wl': ('Same, n=160 with recoil background from '
                                            'E_nu > 2 MeV. From 1_Main.nb.', regions['n160b']),
    }

    testing = header('16_testing', thr, 'Consistency tests. From 4_Testing.nb.') \
        + body_of(P('4_Testing.wl'))

    checks = ''.join(regions['checks']).rstrip()
    checks = (checks + '\n') if checks else '(* none *)\n'
    main = header('00_main', thr,
                  'Orchestrator: loads the numbered modules in order. '
                  'From 1_Main.nb (Setting/Calculations).') \
        + MAIN_TMPL.format(checks=checks)

    # ---------- write everything ----------
    out = {'00_main.wl': main, '01_setup.wl': setup, '02_functions.wl': functions,
           '03_constants_nuflux.wl': constants, '04_exprate_defs.wl': exprate,
           '05_input_data.wl': input_data, '16_testing.wl': testing}
    for fname, (prov, key) in calc_files.items():
        out[fname] = header(fname[:-3], thr, prov) + '\n\n'.join(calc_targets[key]) + '\n'
    for fname, (prov, text) in main_files.items():
        out[fname] = header(fname[:-3], thr, prov) + text.strip() + '\n'

    for fname, text in out.items():
        (src / fname).write_text(text)
        print(f'  wrote {thr}/src/{fname}  ({len(text)} bytes)')

    for old in ['1_Main.wl', '2_1_Input_general.wl', '2_2_Input_special.wl',
                '2_3_Input_Data.wl', '4_Testing.wl', '5_Calculation.wl',
                'Units.wl', 'Functions.wl', 'Constants_NuFlux.wl', 'ExpRate.wl']:
        p = P(old)
        if p.exists():
            p.unlink()
    print(f'  removed legacy-named files in {thr}/src')


for thr in THRESHOLDS:
    print(f'== {thr} ==')
    process(thr)
print('RESTRUCTURE_DONE')
