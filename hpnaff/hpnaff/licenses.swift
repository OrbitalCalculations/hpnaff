//
//  licenses.swift
//  hpnaff
//
//  Created by Heiko Pälike on 08/06/2018.
//  Copyright © 2018 Heiko Pälike. All rights reserved.
//

import Foundation

let licenses = """
The software “hpNAFF” uses libraries and code from the following
projects (and licenses):


1)	SDDS Toolkit
	(https://github.com/epicsdeb/sdds/blob/master/extensions/src/SDDS/LICENSE)
2)	CommandLineKit (https://github.com/jatoben/CommandLine)
3)	Surge: MIT License
	(https://github.com/alejandro-isaza/Upsurge/blob/master/LICENSE)
4)	Rainbow: MIT License
	(https://github.com/onevcat/Rainbow/blob/master/LICENSE)

========================================================================
SDDS Toolkit License: Copyright
(c) 2002 University of Chicago. All rights reserved.

SDDS ToolKit is distributed subject to the following license conditions:

 SOFTWARE LICENSE AGREEMENT Software: SDDS ToolKit

 1. The "Software", below, refers to SDDS ToolKit (in either source
 code, or binary form and accompanying documentation). Each licensee is
 addressed as "you" or "Licensee."

 2. The copyright holders shown above and their third-party licensors
 hereby grant Licensee a royalty-free nonexclusive license, subject to
 the limitations stated herein and U.S. Government license rights.

 3. You may modify and make a copy or copies of the Software for use
 within your organization, if you meet the following conditions: a.
 Copies in source code must include the copyright notice and this
 Software License Agreement. b. Copies in binary form must include the
 copyright notice and this Software License Agreement in the
 documentation and/or other materials provided with the copy.

 4. You may modify a copy or copies of the Software or any portion of
 it, thus forming a work based on the Software, and distribute copies of
 such work outside your organization, if you meet all of the following
 conditions: a. Copies in source code must include the copyright notice
 and this Software License Agreement; b. Copies in binary form must
 include the copyright notice and this Software License Agreement in the
 documentation and/or other materials provided with the copy; c.
 Modified copies and works based on the Software must carry prominent
 notices stating that you changed specified portions of the Software.

 5. Portions of the Software resulted from work developed under a U.S.
 Government contract and are subject to the following license: the
 Government is granted for itself and others acting on its behalf a
 paid-up, nonexclusive, irrevocable worldwide license in this computer
 software to reproduce, prepare derivative works, and perform publicly
 and display publicly.

 6. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT
 WARRANTY OF ANY KIND. THE COPYRIGHT HOLDERS, THEIR THIRD PARTY
 LICENSORS, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY,
 AND THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED,
 INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE, TITLE OR NON-INFRINGEMENT, (2) DO NOT
 ASSUME ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 COMPLETENESS, OR USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT
 USE OF THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4) DO
 NOT WARRANT THAT THE SOFTWARE WILL FUNCTION UNINTERRUPTED, THAT IT IS
 ERROR-FREE OR THAT ANY ERRORS WILL BE CORRECTED.

 7. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT HOLDERS,
 THEIR THIRD PARTY LICENSORS, THE UNITED STATES, THE UNITED STATES
 DEPARTMENT OF ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF ANY KIND OR
 NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF PROFITS OR LOSS OF DATA,
 FOR ANY REASON WHATSOEVER, WHETHER SUCH LIABILITY IS ASSERTED ON THE
 BASIS OF CONTRACT, TORT (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR
 OTHERWISE, EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 POSSIBILITY OF SUCH LOSS OR DAMAGES.


========================================================================
Surge (https://github.com/alejandro-isaza/Upsurge/blob/master/LICENSE)

The MIT License (MIT)

Copyright © 2014-2018 the Surge contributors

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


========================================================================
CommandLineLit (https://github.com/jatoben/CommandLine) License

Copyright (c) 2014 Ben Gollmer.

Licensed under the Apache License, Version 2.0 (the "License"); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.



========================================================================

Rainbow (https://github.com/onevcat/Rainbow/blob/master/LICENSE)

The MIT License (MIT)

Copyright (c) 2015 Wei Wang

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""

let references = """
This is a Swift implementation of https://ops.aps.anl.gov/manuals/SDDStoolkit/SDDStoolkitsu61.html
Related efforts:
https://github.com/MichaelEhrlichman/FortNAFF
https://github.com/adrn/SuperFreq
https://github.com/nkarast/PyNAFF
References:
1. Laskar, J., 1990, The chaotic motion of the Solar System. A numerical estimate of the size of the chaotic zones, Icarus, 88, 266-291.
2. Laskar, J., 1993, Frequency analysis for multi-dimensional systems. Global dynamics and diffusion, Physica D, 67, 257-281.
3. Dumas, S., Laskar, J., 1993, Global Dynamics and Long-Time Stability in Hamiltonian Systems Via Numerical Frequency Analysis, Phys. Rev. Letters, 70 (20), 2975-2979.
4. Laskar, J. : 1999, Introduction to frequency map analysis, in proc. of NATO ASI 533 3DHAM95, S'Agaro, Spain, 134150.
5. Papaphilippou, Y., Frequency maps for LHC models, PAC99.
6. Papahilippou, Y. Zimmermann, F., Weak-strong beam-beam simulations for the Large Hadron Collider, Phys. Rev. ST Accel. Beams 2, 104001 (1999).
7. Robin, D., Steir, C., Laskar, J., Nadolski, L. : 2000, Global dynamics of the ALS revealed through experimental Frequency Map Analysis, Phys. Rev. Let., 85, pp. 558-561.
8. Laskar, J., Frequency map analysis and quasiperiodic decompositions. preprint (https://arxiv.org/pdf/math/0305364.pdf) (2003).
9. Valluri & Merritt, 1998 (http://iopscience.iop.org/article/10.1086/306269/fulltext/37764.text.html#sc2)
"""
