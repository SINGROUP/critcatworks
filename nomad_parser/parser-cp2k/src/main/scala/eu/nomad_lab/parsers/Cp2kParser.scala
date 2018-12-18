/*
 * Copyright 2015-2018 Lauri Himanen, Fawzi Mohamed, Ankit Kariryaa
 * 
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object Cp2kParser extends SimpleExternalParserGenerator(
  name = "Cp2kParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("Cp2kParser")) ::
      ("parserId" -> jn.JString("Cp2kParser" + lab.Cp2kVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.Cp2kVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """  \*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s(?<cp2kStartedAt>.*)
 \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*
 \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*
 \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*
  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*
(?:\s*\n|                                      \s+.*
)*
(?:\s*CP2K\| version string:\s*(?<cp2kVersionString>.*)
)?(?:\s*CP2K\| source code revision number:\s*(?<cp2kRevision>.*)
)?""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/cp2k/parser/parser-cp2k/cp2kparser/scalainterface.py",
    "${mainFilePath}"),
  cmdCwd = "${mainFilePath}/..",
  resList = Seq(
    "parser-cp2k/cp2kparser/__init__.py",
    "parser-cp2k/cp2kparser/setup_paths.py",
    "parser-cp2k/cp2kparser/parser.py",
    "parser-cp2k/cp2kparser/generic/__init__.py",
    "parser-cp2k/cp2kparser/generic/inputparsing.py",
    "parser-cp2k/cp2kparser/versions/__init__.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/__init__.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/singlepointparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/geooptparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/mdparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/singlepointforceparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/inputparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/commonparser.py",
    "parser-cp2k/cp2kparser/versions/cp2k262/input_data/cp2k_input_tree.pickle",
    "parser-cp2k/cp2kparser/scalainterface.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/cp2k.nomadmetainfo.json",
    "nomad_meta_info/cp2k.general.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-cp2k" -> "parsers/cp2k/parser/parser-cp2k",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
