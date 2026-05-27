# Template Switch Aligner

<table border="0">
 <tr>
    <td>Binary:</td>
    <td><a href="https://crates.io/crates/tsalign"><img src="https://img.shields.io/crates/v/tsalign.svg?style=flat-square" alt="Binary version" /></a></td>
    <td><a href="https://crates.io/crates/tsalign"><img src="https://img.shields.io/crates/d/tsalign.svg?style=flat-square" alt="Binary downloads" /></a></td>
    <td></td>
 </tr>
 <tr>
    <td>Library:</td>
    <td><a href="https://crates.io/crates/lib_tsalign"><img src="https://img.shields.io/crates/v/lib_tsalign.svg?style=flat-square" alt="Library version" /></a></td>
    <td><a href="https://crates.io/crates/lib_tsalign"><img src="https://img.shields.io/crates/d/lib_tsalign.svg?style=flat-square" alt="Library downloads" /></a></td>
    <td><a href="https://docs.rs/lib_tsalign"><img src="https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square" alt="Library docs" /></a></td>
 </tr>
 <tr>
    <td>Python bindings:</td>
    <td><a href="https://pypi.org/project/tsalign/"><img src="https://img.shields.io/pypi/v/tsalign?style=flat-square" alt="Python bindings version" /></a></td>
    <td><a href="https://pypi.org/project/tsalign/"><img src="https://img.shields.io/pypi/dm/tsalign?style=flat-square" alt="Python bindings downloads" /></a></td>
    <td></td>
 </tr>
</table>

## Features

* Align two genomic sequences while allowing for template switch mutations (TSMs) based on the [four-point model](https://doi.org/10.1101/gr.214973.116).
* Visualise the alignment with template switching jumps.

## Installation

### Via Cargo (Preferred)

1. Install the rust toolchain by going to [rustup.rs](https://rustup.rs/) and following the instructions.
   Don't worry, on unix-like systems there is just a single command to execute.

2. Run `cargo install tsalign`.

3. You can now run `tsalign` on your command line from anywhere.

If you ever want to update to a new release, simply run `cargo install tsalign` again.

### Python bindings

1. Install the Python package manager [`pip`](https://pypi.org/project/pip/).

2. Run `pip install tsalign`.

3. Refer to [the Python bindings README](python_bindings/README.md) for more information.

### From Source (For Developers)

1. Install the rust toolchain by going to [rustup.rs](https://rustup.rs/) and following the instructions.
   Don't worry, on unix-like systems there will just be a single command to be executed.

2. Clone this git repository using `git clone <url of this repository>`.

3. From within the root of the git repository, you can run `cargo run --release -- align <arguments>` to run tsalign, where `<arguments>` are the arguments that are passed to tsalign.

Cargo acts as a wrapper here, ensuring that whenever you make changes to the code, it will be recompiled if necessary.
Hence, for updating, it is enough to do a `git pull`.

## Usage

For a brief overview of the available subcommands, run `tsalign --help`.
For a brief overview of the available options of each subcommand, run `tsalign <subcommand> --help`.

### Aligning

Tsalign takes as input either a single fasta file with two fasta records, or two separate fasta files.
The output is a toml file which contains alignment statistics and can be used to visualise the alignment.

The alignment metric can be configured via a config file named `config.tsa`.
The path to the file can be passed via the `-c` parameter.
Only pass the directory that contains the file, and not the file itself.

Note that tsalign is very resource-intensive, so it often makes sense to use the `--memory-limit` parameter, which takes a memory limit in bytes.
This limit is applied approximately, with an accuracy of about a factor of two.

If you have two sequences in a file `pair.fa`, computing an alignment could look as follows.

```bash
tsalign align -p pair.fa -o alignment.toml -c sample_tsa_config --memory-limit 1000000000
```

#### Choosing a different alphabet

If you want to align strings with an alphabet other than the default DNA alphabet, consider using the parameter `-a`.
Remember to also update the `config.tsa` file with a metric that corresponds to the characters available in the alphabet.
Available alphabets are:

* **dna**: ACGT
* **dna-n**: ACGTN
* **rna**: ACGU
* **rna-n**: ACGUN
* **dna-iupac**: ABCDGHKMNRSTVWY
* **rna-iupac**: ABCDGHKMNRSUVWY

#### Skip characters

If your input file(s) contain characters outside of the alphabet you specified, such as `-` characters marking an existing alignment, you must mark them as ignored using the `--skip-characters` parameter.
This does not pertain the special `|` character if `--use-embedded-rq-ranges` is used as explained below.

#### Disabling template switches

If you want to compute an alignment without template switches, you can use the `--no-ts` parameter.
This can be useful for comparing the optimal template switching alignment with the optimal alignment without template switches.

### Modifying the cost function

We designed tsalign to allow for very detailed alignment costs regarding not just the gap-affine alignment, but also the geometry of a TSM.
Consult [`sample_tsa_config/config.tsa`](sample_tsa_config/config.tsa) for an example.
The file format is very strict in that the order of options must stay the same.
We explain the content of the file step-by-step.

#### Using different gap-affine costs for TSM flanks

You can use different gap-affine costs for the upstream and downstream sequences (flanks) of TSMs.
Increasing the size of these flanks affects performance a lot.
The flank lengths are set here, and the different costs are set at the end of the file.

```txt
# Limits

left_flank_length = 0
right_flank_length = 0
```

#### TSM base cost

Each TSM incurs a cost based on its type.
The type is written as `<descendant><ancestor><direction>`.
In descendant and ancestor, the letter `r` refers to the first sequence (reference), and the letter `q` refers to the second sequence (query).
In direction, the letter `f` refers to a repeat, and the letter `r` refers to a TSM.
Repeats are not well supported, so consider giving them cost `inf` to disable them entirely.

```txt
# Base Cost

rrf_cost = inf
rqf_cost = inf
qrf_cost = inf
qqf_cost = inf
rrr_cost = 3
rqr_cost = 2
qrr_cost = 2
qqr_cost = 3
```

#### TSM jump costs

Each TSM additionally incurs cost based on its geometry.
The costs are a piecewise constant function, where the first row is the first input value that the constant cost applies to, and the second row is the constant cost.
`RQQROffset`, `RRQQOffset`, and `LengthDifference` must be 0-centered V-shaped, i.e. non-ascending for inputs `< 0` and non-descending for inputs `>= 0`.

`RQQROffset` and `RRQQOffset` are costs based on the length of the 1-2-jump of the TSM.
`RQQROffset` is applied to TSMs where ancestor and descendant are different, while `RRQQOffset` is applied to TSMs where ancestor and descendant are the same.
`Length` is the cost based on the length of the 2-3-alignment of the TSM.
`LengthDifference` is the cost based on the difference between the length of the 2-3-alignment and the difference between points 1 and 4.
`ForwardAntiDescendantGap` is the cost based on the difference between points 1 and 4, specifically `SP4 - SP1`.

```txt
# Jump Costs

RQQROffset
 -inf -100 101
  inf    0 inf

RRQQOffset
 -inf -100 101
  inf    0 inf

Length
   0 5 6 7 8 100
 inf 5 3 1 0 inf

LengthDifference
 -inf -100 101
  inf    0 inf

ForwardAntiDescendantGap
 -inf   1
    0 inf

ReverseAntiDescendantGap
 -inf
    0
```

#### Gap-affine edit costs

Different gap-affine edit costs are applied to different parts of the alignment.
Outside of any TSM, the primary costs are applied.
The `SubstitutionCostTable`, `GapOpenCostVector`, and `GapExtendCostVector` must have the letters of the chosen alphabet.
The example is made for `dna-n`.
For e.g. `dna`, the rows and columns with `N` must be removed, while for other alphabets, other columns need to be removed or added.
The costs are gap-affine, where the first character of a gap is priced with `GapOpenCostVector`, and all following characters are priced with `GapExtendCostVector`.

```txt
# Primary Edit Costs

SubstitutionCostTable
  |  A  C  G  T  N
--+---------------
A |  0  2  2  2  0
C |  2  0  2  2  0
G |  2  2  0  2  0
T |  2  2  2  0  0
N |  0  0  0  0  0

GapOpenCostVector
 A C G T N
 3 3 3 3 3

GapExtendCostVector
 A C G T N
 1 1 1 1 1
```

Following the primary edit costs, there are analogous sections.
`Secondary Forward Edit Costs` are the costs for the 2-3-alignment of a repeat.
`Secondary Reverse Edit Costs` are the costs for the 2-3-alignment of a TSM.
`Left Flank Edit Costs` are the costs in the upstream flank of a TSM.
`Right Flank Edit Costs` are the costs in the downstream flank of a TSM.

```txt
# Secondary Forward Edit Costs

# Secondary Reverse Edit Costs

# Left Flank Edit Costs

# Right Flank Edit Costs
```

### SVG/PNG Visualisation

To visualise the alignment, you can use the `tsalign show` subcommand.
For an optimal experience, compute both the alignment with template switches and the alignment without template switches.
For example:

```bash
tsalign align -p pair.fa -o alignment.toml -c sample_tsa_config --memory-limit 1000000000
tsalign align -p pair.fa -o alignment-no-ts.toml -c sample_tsa_config --memory-limit 1000000000 --no-ts
```

Then, create the visualisation in SVG format as follows:

```bash
tsalign show -i alignment.toml -n alignment-no-ts.toml -s alignment.svg
```

Since SVGs are not always well supported, you can also use the switch `-p` to render the visualisation also as PNG.
Note that the `-p` switch requires setting the `-s` parameter.

```bash
tsalign show -i alignment.toml -n alignment-no-ts.toml -p -s visualisation.svg
```

There are additional switches and parameters to influence the visualisation.

* `-a` adds arrows for the jumps of the alignments
* `-c` renders more of the complements of the input sequences
* `-e` renders the heuristically computed uncertainties of a TSM with grey characters
* `-z X` renders only `X` bases of context around the TSMs

```bash
tsalign show -i alignment.toml -n alignment-no-ts.toml -ps visualisation.svg -ace -z 5
```

### Setting the alignment range

When searching for template switches, researchers are often interested in specific mutation hotspots that are only tens of characters long, while the 2-3-alignment of the template switch may also align outside of the mutation hotspot.
This requires to pass longer sequences to tsalign, to make the context around the mutation hotspot available for alignment.
However, since tsalign's performance is very sensitive to the lengths of the input sequences, it is important to specify the range of the mutation hotspot, such that tsalign computes the alignment only of this shorter region, while only allowing 2-3-alignments of template switches to align outside of the region.

For an example, consider the following input sequences:

Original `input.fa`.

```fasta
>reference
ACACACCCAACGCGGG
>query
ACAAACGTGTCGCGCG
```

We want to check if the substitution of `CCAA` to `GTGT` can be better explained with a template switch.
Therefore, we want tsalign to focus on this region.
We can insert `|` characters to mark the focus region, including an extra character at the beginning and the end for good measure.

Modified `input.delimited.fa`.

```fasta
>reference
ACACA|CCCAAC|GCGGG
>query
ACAAA|CGTGTC|GCGCG
```

Now we can run tsalign as follows:

```bash
tsalign align -p input.delimited.fa --use-embedded-rq-ranges
```

Now, tsalign will use much less resources, as it can ignore the non-matches before and after the focus region.
