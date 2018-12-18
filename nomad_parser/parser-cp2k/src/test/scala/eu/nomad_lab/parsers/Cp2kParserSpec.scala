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

import org.specs2.mutable.Specification

object Cp2kParserSpec extends Specification {
  "Cp2kParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(Cp2kParser, "parsers/cp2k/test/examples/energy_force/si_bulk8.out", "json-events") must_== ParseResult.ParseSuccess
    }
  }

  "test energy_force with json" >> {
    ParserRun.parse(Cp2kParser, "parsers/cp2k/test/examples/energy_force/si_bulk8.out", "json") must_== ParseResult.ParseSuccess
  }
  "test geo_opt with json" >> {
    ParserRun.parse(Cp2kParser, "parsers/cp2k/test/examples/geo_opt/H2O.out", "json") must_== ParseResult.ParseSuccess
  }
  "test md with json" >> {
    ParserRun.parse(Cp2kParser, "parsers/cp2k/test/examples/md/H2O-32.out", "json") must_== ParseResult.ParseSuccess
  }
}
