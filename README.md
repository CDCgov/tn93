# Tamura-Nei Distance Calculation with Python

## Overview

This repository contains a Python implementation of the Tamura-Nei (TN93) distance calculation. When provided with two sequences to compare, this software returns a distance value between 0 and 1 using the Tamura-Nei nucleotide substitution model. See  [Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees](https://pubmed.ncbi.nlm.nih.gov/8336541/) and [NOTES](https://github.com/CDCgov/tn93/blob/main/NOTES.md) for more information on the algorithm.

## Usage

This tool is primarily meant as a library to be imported and used in custom analysis code, but can also be used to directly calculate the pairwise distances for a set of sequences in a FASTA file.

First, install using pip

```bash
pip install tn93
```

or clone this respository and copy `src/tn93/tn93.py` to your working directory. To calculate the distance between a pair of sequences,

```python
from Bio import SeqIO
import tn93
# Read in a FASTA file to get sequences
seqs = [ x for x in SeqIo.parse("your_sequences.fasta", format="fasta") ]
tn93 = tn93.TN93()
distance = tn93.tn93_distance(seqs[0], seqs[1], "RESOLVE")
```

Alternatively, the module can be run from the command line and provided with a sequence file and match mode to produce a JSON file with the pairwise distances.

```bash
python tn93.py --input_file example_seqs.fasta --match_mode RESOLVE --output example_seqs_resolve_distance.json
```

There are four distinct match modes:

* SKIP, which ignores ambiguous positions
* GAPMM, which treats gaps appearing in only one sequence as mismatches
* AVERAGE, which takes the average of the possible resolution values
* RESOLVE, which tries to resolve the ambiguity to a single nucleotide, averages if that fails

## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
