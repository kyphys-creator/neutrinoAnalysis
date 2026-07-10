#!/usr/bin/env wolframscript

(* Convert a single .nb file to .wl by extracting Code/Input cell contents.

   String cells are emitted verbatim. BoxData cells are converted to InputForm
   text via MakeExpression + HoldForm.

   Usage:  wolframscript -file convert_nb_to_wl.wls <input.nb> <output.wl>

   Replaces NotebookDirectory[] with DirectoryName[$InputFileName] in the
   Setup cell. Updates Get[FileNameJoin[{nbDir, "XX.nb"}]] → ...wl path.
*)

args = Rest[$ScriptCommandLine];
If[Length[args] < 2,
  Print["usage: convert_nb_to_wl.wls <input.nb> <output.wl>"];
  Exit[1]
];

inputPath  = args[[1]];
outputPath = args[[2]];

stripHoldForm[s_String] := StringReplace[
  s,
  {StartOfString ~~ "HoldForm[" -> "", "]" ~~ EndOfString -> ""}
];

(* Convert one box statement to InputForm text *)
$failedStmts = 0;
boxStmtToText[box_] := Module[{expr, text, clean},
  (* raw "\n" tokens inside RowBox lists break MakeExpression; treat them as
     plain indenting newlines *)
  clean = box /. "\n" -> "\[IndentingNewLine]";
  expr = Quiet@MakeExpression[clean, StandardForm];
  (* some box constructs (nested compound definitions) fail MakeExpression but
     parse fine via ToExpression with a HoldComplete wrapper *)
  If[Head[expr] =!= HoldComplete,
    expr = Quiet@ToExpression[clean, StandardForm, HoldComplete]];
  If[Head[expr] =!= HoldComplete, $failedStmts++; Return[""]];
  text = ToString[HoldForm @@ expr, InputForm, PageWidth -> 120];
  StringTrim @ stripHoldForm[text]
];

(* Convert a cell's content to text *)
cellToText[s_String] := s;
cellToText[BoxData[boxes_]] := Module[{stmts, results},
  stmts = If[Head[boxes] === List,
    DeleteCases[boxes, _String?(StringMatchQ[#, Whitespace | "\n"] &)],
    {boxes}
  ];
  results = Select[boxStmtToText /@ stmts, # =!= "" && # =!= "Null" &];
  StringRiffle[results, "\n"]
];
cellToText[_] := "";  (* fallback *)

(* Update Setup pattern: NotebookDirectory[] → DirectoryName[$InputFileName] *)
(* Also rewrite Get[FileNameJoin[{nbDir, "XX.nb"}]] → ...wl              *)
postProcess[code_String] := StringReplace[code,
  {
    (* MakeExpression artifacts: cells with embedded comments produce
       (EXPR; )*Null  or  Null*((s1; )*Null*(s2; ))  — clean these up.
       The product-splitting rule REQUIRES the "; )" ending so that genuine
       multiplications like (a - b)*(c*d) are never touched; only compound-
       statement products (s1; )*(s2; ) — the residue of stripped comments —
       are split into separate statements. *)
    RegularExpression["\\(([^()]+?;)\\s*\\)(\\*Null)+"] :> "$1",
    RegularExpression["Null\\*\\(([^()]+?;)\\s*\\)"]    :> "$1",
    RegularExpression["(;\\s*)\\)\\*\\("]               :> "$1); (",
    "*Null"                                              -> "",
    "Null*"                                              -> "",

    (* Setup pattern *)
    "NotebookDirectory[]"            -> "DirectoryName[$InputFileName]",
    "nbDir = DirectoryName[$InputFileName]" -> "fileDir = DirectoryName[$InputFileName]",
    RegularExpression["nbDir"]       -> "fileDir",

    (* Get-chain *)
    RegularExpression["\"(\\d\\d_[A-Za-z_]+)\\.nb\""] -> "\"$1.wl\""
  }
];

Print["Reading ", inputPath, "..."];
nb = Get[inputPath];
If[Head[nb] =!= Notebook,
  Print["Error: not a notebook expression"];
  Exit[1]
];

(* Walk cells in order, keeping section headings as comment markers so the
   output can be split by topic afterwards. *)
headingText[s_String] := s;
headingText[TextData[l_List]] := StringJoin[Cases[l, _String, Infinity]];
headingText[TextData[s_String]] := s;
headingText[_] := "";

cells = Cases[nb,
  Cell[content_, style : ("Title" | "Section" | "Subsection" | "Subsubsection" |
       "Code" | "Input" | "InitialFunc" | "InitializeFunction"), ___] :>
    {style, content},
  Infinity];
Print["  ", Length[cells], " cells (incl. headings)"];

toText[{style_, content_}] := Which[
  MemberQ[{"Title", "Section", "Subsection", "Subsubsection"}, style],
  Module[{h = StringTrim@StringReplace[headingText[content], "\n" -> " "]},
    If[h === "", "",
      "(* ===== " <> style <> ": " <> h <> " ===== *)"]],
  True, cellToText[content]
];

texts = toText /@ cells;
combined = StringRiffle[Select[texts, # =!= "" &], "\n\n"];
combined = postProcess[combined];

(* Add a header comment *)
header = StringTemplate[
  "(* ========================================================================== *)
(*  `name`\n*)
(*  Auto-converted from `srcname` by convert_nb_to_wl.wls\n*)
(* ========================================================================== *)\n\n"
][<|"name" -> FileBaseName[outputPath], "srcname" -> FileNameTake[inputPath]|>];

Export[outputPath, header <> combined, "Text"];
Print["  failed statements: ", $failedStmts];
Print["Wrote ", outputPath, " (", FileByteCount[outputPath], " bytes)"];
