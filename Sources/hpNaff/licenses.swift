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
5)  a modified version of Polynomial from
  damuellen/Utilities/Physics/Coefficients+Ratio.swift: Apache License v2
  (https://github.com/damuellen/Utilities/blob/main/Sources/Physics/Coefficients%2BRatio.swift)
  Modification: fitting routine from https://github.com/natedomin/polyfit
6)  Doggie: MIT License
  (https://github.com/SusanDoggie/Doggie/blob/main/LICENSE)
7)  Swift-numerics: Apache License v2
(https://github.com/apple/swift-numerics/blob/main/LICENSE.txt)


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

========================================================================

damuellen/Utilities/Physics/Coefficients+Ratio.swift: Apache License v2
(https://github.com/damuellen/Utilities/blob/main/Sources/Physics/Coefficients%2BRatio.swift)

Copyright 2017 Daniel Müllenborn

Apache License
Version 2.0, January 2004
http://www.apache.org/licenses/

TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

1. Definitions.

"License" shall mean the terms and conditions for use, reproduction, and distribution as defined by Sections 1 through 9 of this document.

"Licensor" shall mean the copyright owner or entity authorized by the copyright owner that is granting the License.

"Legal Entity" shall mean the union of the acting entity and all other entities that control, are controlled by, or are under common control with that entity. For the purposes of this definition, "control" means (i) the power, direct or indirect, to cause the direction or management of such entity, whether by contract or otherwise, or (ii) ownership of fifty percent (50%) or more of the outstanding shares, or (iii) beneficial ownership of such entity.

"You" (or "Your") shall mean an individual or Legal Entity exercising permissions granted by this License.

"Source" form shall mean the preferred form for making modifications, including but not limited to software source code, documentation source, and configuration files.

"Object" form shall mean any form resulting from mechanical transformation or translation of a Source form, including but not limited to compiled object code, generated documentation, and conversions to other media types.

"Work" shall mean the work of authorship, whether in Source or Object form, made available under the License, as indicated by a copyright notice that is included in or attached to the work (an example is provided in the Appendix below).

"Derivative Works" shall mean any work, whether in Source or Object form, that is based on (or derived from) the Work and for which the editorial revisions, annotations, elaborations, or other modifications represent, as a whole, an original work of authorship. For the purposes of this License, Derivative Works shall not include works that remain separable from, or merely link (or bind by name) to the interfaces of, the Work and Derivative Works thereof.

"Contribution" shall mean any work of authorship, including the original version of the Work and any modifications or additions to that Work or Derivative Works thereof, that is intentionally submitted to Licensor for inclusion in the Work by the copyright owner or by an individual or Legal Entity authorized to submit on behalf of the copyright owner. For the purposes of this definition, "submitted" means any form of electronic, verbal, or written communication sent to the Licensor or its representatives, including but not limited to communication on electronic mailing lists, source code control systems, and issue tracking systems that are managed by, or on behalf of, the Licensor for the purpose of discussing and improving the Work, but excluding communication that is conspicuously marked or otherwise designated in writing by the copyright owner as "Not a Contribution."

"Contributor" shall mean Licensor and any individual or Legal Entity on behalf of whom a Contribution has been received by Licensor and subsequently incorporated within the Work.

2. Grant of Copyright License. Subject to the terms and conditions of this License, each Contributor hereby grants to You a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable copyright license to reproduce, prepare Derivative Works of, publicly display, publicly perform, sublicense, and distribute the Work and such Derivative Works in Source or Object form.

3. Grant of Patent License. Subject to the terms and conditions of this License, each Contributor hereby grants to You a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable (except as stated in this section) patent license to make, have made, use, offer to sell, sell, import, and otherwise transfer the Work, where such license applies only to those patent claims licensable by such Contributor that are necessarily infringed by their Contribution(s) alone or by combination of their Contribution(s) with the Work to which such Contribution(s) was submitted. If You institute patent litigation against any entity (including a cross-claim or counterclaim in a lawsuit) alleging that the Work or a Contribution incorporated within the Work constitutes direct or contributory patent infringement, then any patent licenses granted to You under this License for that Work shall terminate as of the date such litigation is filed.

4. Redistribution. You may reproduce and distribute copies of the Work or Derivative Works thereof in any medium, with or without modifications, and in Source or Object form, provided that You meet the following conditions:

You must give any other recipients of the Work or Derivative Works a copy of this License; and
You must cause any modified files to carry prominent notices stating that You changed the files; and
You must retain, in the Source form of any Derivative Works that You distribute, all copyright, patent, trademark, and attribution notices from the Source form of the Work, excluding those notices that do not pertain to any part of the Derivative Works; and
If the Work includes a "NOTICE" text file as part of its distribution, then any Derivative Works that You distribute must include a readable copy of the attribution notices contained within such NOTICE file, excluding those notices that do not pertain to any part of the Derivative Works, in at least one of the following places: within a NOTICE text file distributed as part of the Derivative Works; within the Source form or documentation, if provided along with the Derivative Works; or, within a display generated by the Derivative Works, if and wherever such third-party notices normally appear. The contents of the NOTICE file are for informational purposes only and do not modify the License. You may add Your own attribution notices within Derivative Works that You distribute, alongside or as an addendum to the NOTICE text from the Work, provided that such additional attribution notices cannot be construed as modifying the License.

You may add Your own copyright statement to Your modifications and may provide additional or different license terms and conditions for use, reproduction, or distribution of Your modifications, or for any such Derivative Works as a whole, provided Your use, reproduction, and distribution of the Work otherwise complies with the conditions stated in this License.
5. Submission of Contributions. Unless You explicitly state otherwise, any Contribution intentionally submitted for inclusion in the Work by You to the Licensor shall be under the terms and conditions of this License, without any additional terms or conditions. Notwithstanding the above, nothing herein shall supersede or modify the terms of any separate license agreement you may have executed with Licensor regarding such Contributions.

6. Trademarks. This License does not grant permission to use the trade names, trademarks, service marks, or product names of the Licensor, except as required for reasonable and customary use in describing the origin of the Work and reproducing the content of the NOTICE file.

7. Disclaimer of Warranty. Unless required by applicable law or agreed to in writing, Licensor provides the Work (and each Contributor provides its Contributions) on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE. You are solely responsible for determining the appropriateness of using or redistributing the Work and assume any risks associated with Your exercise of permissions under this License.

8. Limitation of Liability. In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.

9. Accepting Warranty or Additional Liability. While redistributing the Work or Derivative Works thereof, You may choose to offer, and charge a fee for, acceptance of support, warranty, indemnity, or other liability obligations and/or rights consistent with this License. However, in accepting such obligations, You may act only on Your own behalf and on Your sole responsibility, not on behalf of any other Contributor, and only if You agree to indemnify, defend, and hold each Contributor harmless for any liability incurred by, or claims asserted against, such Contributor by reason of your accepting any such warranty or additional liability.

END OF TERMS AND CONDITIONS

========================================================================
Doggie: (https://github.com/SusanDoggie/Doggie)

The MIT License (MIT)

Copyright (c) 2015 - 2020 SusanDoggie. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


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
