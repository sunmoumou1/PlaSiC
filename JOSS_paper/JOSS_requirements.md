# Submitting a paper to JOSS

If you've already developed a fully featured research code, released it under an [OSI-approved license](https://opensource.org/licenses), and written good documentation and tests, then we expect that it should take perhaps an hour or two to prepare and submit your paper to JOSS.
But please read these instructions carefully for a streamlined submission.

## Submission requirements

- The software must be open source as per the [OSI definition](https://opensource.org/osd).
- The software must be hosted at a location where users can browse the source code files, open issues, and propose code changes without manual approval of (or payment for) accounts
- The software must have an **obvious** research application.
- You must be a major contributor to the software you are submitting, and have a GitHub account to participate in the review process.
- Your paper must not focus on new research results accomplished with the software.
- Your paper (`paper.md` and BibTeX files, plus any figures) must be hosted in a Git-based repository together with your software.
- The paper may be in a short-lived branch which is never merged with the default, although if you do this, make sure this branch is _created_ from the default so that it also includes the source code of your submission.

In addition, the software associated with your submission must:

- Be stored in a repository that can be cloned without registration.
- Be stored in a repository where the software source files are browsable online without registration.
- Have an issue tracker that is readable without registration.
- Permit individuals to create issues/file tickets against your repository.

### What we mean by research software

JOSS publishes articles about research software. This definition includes software that: solves complex modeling problems in a scientific context (physics, mathematics, biology, medicine, social science, neuroscience, engineering); supports the functioning of research instruments or the execution of research experiments; extracts knowledge from large data sets; offers a mathematical library, or similar. While useful for many areas of research, pre-trained machine learning models and notebooks are not in-scope for JOSS. 

### Scope and significance

JOSS publishes articles about software that has demonstrated clear research impact. Your software should represent a meaningful contribution to the research community rather than being a one-off tool for a single analysis. Whether or not you use AI assistance in your development, submissions should demonstrate the irreplaceable human contributions: problem framing, key design decisions, thoughtful architectural choices, and practices that make software usable and sustainable for others.

Some factors that may be considered by editors and reviewers when evaluating scope and significance include:

- **Research impact:** Evidence of publications or analyses using the software, external adopters or integrations.
- **Design thinking:** Meaningful architectural decisions and trade-offs considered. We particularly value work that builds upon or extends existing software ecosystems rather than reinventing solutions where quality alternatives already exist.
- **Open development practices:** Sustained development over time with evidence of collaborative effort, public development history, comprehensive testing, clear documentation, and pathways for community contribution.
- Whether a potential user can easily install, understand, and test the software. (If your software is new, please be sure a colleague has tried it.)
- Whether the software is sufficiently useful that it is _likely to be cited_ by other researchers in your domain.

Projects developed privately are not eligible until there is a public record of open development: at least six months of public history prior to submission, with evidence of releases, public issues and pull requests. A history of contributions and engagement from individuals beyond the original team, across organizations, is especially welcome, though not essential. We particularly value projects that show evidence of open development practices from early stages, rather than private development followed by public release.

In addition, JOSS requires that software should be feature-complete (i.e., no half-baked solutions), packaged appropriately according to common community standards for the programming language being used (e.g., [Python](https://packaging.python.org), [R](https://r-pkgs.org/index.html)), and designed for maintainable extension (not one-off modifications of existing tools). "Minor utility" packages, including "thin" API clients, and single-function packages are not acceptable.

### Pre-review screening criteria

Before entering review, all submissions are evaluated against the following gates. A submission that fails any one of these will receive a desk rejection, but can be resubmitted once the gaps are addressed.

#### Must meet (all required)

**1. Sufficient public development history**

The repository must have been public for more than six months prior to submission, with active development spanning that period. A repository made public immediately before submission, or one showing development concentrated into a few days or weeks,
will not be accepted. We run automated checks on commit distribution — a repo dump is not a history.

**2. Demonstrated research impact**

There **must be evidence that the software is being used for research** — at minimum by the developers themselves, and ideally by others. Acceptable signals include: references in published papers or preprints, documented adoption by other research groups, or clear integration into research workflows. Adoption in currently private workflows (e.g., not shared publicly before an experiment's first papers and/or data releases) is acceptable but must be demonstrated to the editorial team.

Aspirational statements about future use are not sufficient; JOSS will not publish papers that are meant to advertise software that is not yet being used in research.

**3. Good open source practices**

The repository must show active use of open-source workflows. For multi-author projects this means evidence of issues, pull requests, and public discussion. For single-author projects, this bar can be met more broadly — but *multiple* indicators must be present at submission time: a meaningful public commit history over time, tagged releases or a changelog, tests and CI, clear documentation, a CONTRIBUTING file, and stated support or governance expectations. A solo project that is otherwise clearly open and well-maintained will not be rejected solely for lacking a PR workflow. However, a single-author project with *none* of these signals is not ready.

**4. Iterative development over time**

The development history must show ongoing iteration, not a single burst of commits. We look for evidence that the software has been refined through use and feedback over time. A repository where all significant work was added in a concentrated window is a signal that the project was not developed iteratively.

#### Strong positive signal (not a gate, but counts in your favour)

**Community engagement beyond the authors:** non-author issues, pull requests, or discussions are a strong signal of a healthy open project. Submissions with this kind of engagement are well-positioned for review.

#### If your submission doesn't meet these criteria yet

A desk rejection on these grounds may be a "not yet," rather than a "never." If your submission is otherwise solid, we will attempt to write a rejection notice that tells you specifically what is missing. Fix the gaps, continue developing openly, and resubmit in six months or more.

#### AI Usage Policy

**Author use:** The use of generative AI is permitted for most aspects of a JOSS submission (e.g., software creation and review, generating documentation, assisting with paper authoring), however all such use must be disclosed in an "AI usage disclosure" statement which includes:

- **Tool use:** The tools/models used (and versions) and where they were used (code, paper text, docs).
- **The nature and scope of assistance:** e.g., code generation, refactoring, test scaffolding, copy-editing, drafting.
- **Confirmation of review:** Authors must assert that human authors reviewed, edited, validated all AI-assisted outputs and made the core design decisions.

AI is not allowed for conversational interactions between authors and editors or reviewers unless it is being used for translation purposes.

Authors remain fully responsible for the accuracy, originality, licensing, and ethical/legal compliance of all submitted materials. Failure to provide a complete and accurate disclosure of AI usage may be considered an ethical breach. Consequences can include desk rejection, mandatory revisions, and post-publication correction or withdrawal. In cases of intentional misrepresentation or non-disclosure, JOSS reserves the right to notify the authors' institutions, funders, and/or relevant professional or scholarly societies in accordance with standard research-integrity practices.

**Reviewer use:** Reviewers may use generative AI tools to assist with non-substantive tasks (e.g., grammar checks of their review text, contextualizing public materials and code snippets) and must disclose this briefly at the end of their review.

- **Accountability:** The reviewer – not the AI – owns the evaluation. All judgements, recommendations, and technical claims must be the reviewer's and verified by them.
- **Human-only judgements:** When reviewing JOSS submissions (code or papers), all evaluative decisions – scoring, accept/reject recommendations, and assessments of originality, novelty, correctness, significance, and policy/ethics compliance – must be made by the human reviewer. AI tools may assist with analysis, but they must not render or determine any verdicts.


### Co-publication of science, methods, and software

Sometimes authors prepare a JOSS publication alongside a contribution describing a science application, details of algorithm development, and/or methods assessment. In this circumstance, JOSS considers submissions for which the implementation of the software itself reflects a substantial scientific effort. This may be represented by the design of the software, the implementation of the algorithms, creation of tutorials, or any other aspect of the software. We ask that authors indicate whether related publications (published, in review, or nearing submission) exist as part of submitting to JOSS.

#### Other venues for reviewing and publishing software packages

Authors wishing to publish software deemed out of scope for JOSS have a few options available to them:

- Follow [GitHub's guide](https://guides.github.com/activities/citable-code/) on how to create a permanent archive and DOI for your software. This DOI can then be used by others to cite your work.
- Enquire whether your software might be considered by communities such as [rOpenSci](https://ropensci.org) and [pyOpenSci](https://pyopensci.org).

### Should I write my own software or contribute to an existing package?

While we are happy to review submissions in standalone repositories, we also review submissions that are significant contributions made to existing packages. It is often better to have an integrated library or package of methods than a large number of single-method packages.

## Policies

**Disclosure:** All authors must disclose any potential conflicts of interest related to the research in their manuscript, including financial, personal, or professional relationships that may affect their objectivity. This includes any financial relationships, such as employment, consultancies, honoraria, stock ownership, or other financial interests that may be relevant to the research.

**Acknowledgement:** Authors should acknowledge all sources of financial support for the work and include a statement indicating whether or not the sponsor had any involvement in it.

**Review process:** Editors and reviewers must be informed of any potential conflicts of interest before reviewing the manuscript to ensure unbiased evaluation of the research.

**Compliance:** Authors who fail to comply with the COI policy may have their manuscript rejected or retracted if a conflict is discovered after publication.

**Review and Update:** This COI policy will be reviewed and updated regularly to ensure it remains relevant and effective.

### Conflict of Interest policy for authors

An author conflict of interest (COI) arises when an author has financial, personal, or other interests that may influence their research or the interpretation of its results. In order to maintain the integrity of the work published in JOSS, we require that authors disclose any potential conflicts of interest at submission time.

### Preprint Policy

Authors are welcome to submit their papers to a preprint server ([arXiv](https://arxiv.org/), [bioRxiv](https://www.biorxiv.org/), [SocArXiv](https://socopen.org/), [PsyArXiv](https://psyarxiv.com/) etc.) at any point before, during, or after the submission and review process.

Submission to a preprint server is _not_ considered a previous publication.

### Authorship

Purely financial (such as being named on an award) and organizational (such as general supervision of a research group) contributions are not considered sufficient for co-authorship of JOSS submissions, but active project direction and other forms of non-code contributions are. The authors themselves assume responsibility for deciding who should be credited with co-authorship, and co-authors must always agree to be listed. In addition, co-authors agree to be accountable for all aspects of the work, and to notify JOSS if any retraction or correction of mistakes are needed after publication.

### Submissions using proprietary languages/development environments

We strongly prefer software that doesn't rely upon proprietary (paid for) development environments/programming languages. However, provided _your submission meets our requirements_ (including having a valid open source license) then we will consider your submission for review. Should your submission be accepted for review, we may ask you, the submitting author, to help us find reviewers who already have the required development environment installed.

## Submission Process

Before you submit, you should:

- Make your software available in an open repository (GitHub, Bitbucket, etc.) and include an [OSI approved open source license](https://opensource.org/licenses).
- Make sure that the software complies with the [JOSS review criteria](review_criteria). In particular, your software should be full-featured, well-documented, and contain procedures (such as automated tests) for checking correctness.
- Write a short paper in Markdown format using `paper.md` as file name, including a title, summary, author names, affiliations, and key references. See our [example paper](example_paper) to follow the correct format.
- (Optional) create a metadata file describing your software and include it in your repository. We provide [a script](https://gist.github.com/arfon/478b2ed49e11f984d6fb) that automates the generation of this metadata.

### Submitting your paper

Submission is as simple as:

- Filling in the [short submission form](http://joss.theoj.org/papers/new)
- Waiting for the managing editor to start a pre-review issue over in the JOSS reviews repository: https://github.com/openjournals/joss-reviews

### No submission fees

There are no fees for submitting or publishing in JOSS. You can read more about our [cost and sustainability model](http://joss.theoj.org/about#costs).


## Review Process

After submission:

- An Associate Editor-in-Chief will carry out an initial check of your submission, and proceed to assign a handling editor.
- The handling editor will assign two or more JOSS reviewers, and the review will be carried out in the [JOSS reviews repository](https://github.com/openjournals/joss-reviews).
- Authors will respond to reviewer-raised issues (if any are raised) on the submission repository's issue tracker. Reviewer and editor contributions, like any other contributions, should be acknowledged in the repository. 
- **JOSS reviews are iterative and conversational in nature.** Reviewers are encouraged to post comments/questions/suggestions in the review thread as they arise, and authors are expected to respond in a timely fashion.
- **Authors are expected to respond to reviewer feedback promptly.** We generally ask that authors respond to reviewer comments and questions within 2 weeks, and complete requested changes within 4-6 weeks, unless otherwise negotiated with the editor. If you need more time (for example, if the changes requested are particularly substantial), please communicate this to the reviewers and editor on the review thread. Prolonged unresponsiveness may result in the paper being rejected due to lack of engagement.
- Authors and reviewers are asked to be patient when waiting for a response from an editor. Please allow a week for an editor to respond to a question before prompting them for further action.
- Upon successful completion of the review, authors will make a tagged release of the software, and deposit a copy of the repository with a data-archiving service such as [Zenodo](https://zenodo.org/) or [figshare](https://figshare.com/), get a DOI for the archive, and update the review issue thread with the version number and DOI.
- After we assign a DOI for your accepted JOSS paper, its metadata is deposited with CrossRef and listed on the JOSS website.
- The review issue will be closed, and an automatic post from [@JOSS at Mastodon](https://fosstodon.org/@joss) will announce it!


---


# JOSS Paper Format

Submitted articles must use Markdown and must provide a metadata section at the beginning of the article. Format metadata using YAML, a human-friendly data serialization language. The information provided is included in the title and sidebar of the generated PDF. 

---

## What should my paper contain?

```{important}
Begin your paper with a summary of the high-level functionality of your software for a non-specialist reader. Avoid jargon in this section.
```

The paper should be between 750-1750 words. Authors submitting papers significantly longer than 1750 words may be asked to reduce the length of their paper.

Your paper must include the following **required sections**:

- **Summary**: A description of the high-level functionality and purpose of the software for a diverse, *non-specialist audience*.
- **Statement of need**: A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
- **State of the field**: A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear "build vs. contribute" justification explaining your unique scholarly contribution and why existing alternatives are insufficient.
- **Software design**: An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description.
- **Research impact statement**: Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational.
- **AI usage disclosure**: Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified.

Your paper must also include:

- A list of the authors of the software and their affiliations, using the correct format (see the example below).
- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
- Acknowledgement of any financial support.

As this list shows, JOSS papers are expected to contain a limited set of metadata (see example below), and sections for Summary, Statement of need, State of the field, Software design, Research impact statement, AI usage disclosure, Acknowledgements, and References. You can look at an [example accepted paper](example_paper). 

```{important}
Your paper will be reviewed by two or more reviewers in a public GitHub issue. Take a look at the [review checklist](review_checklist) and  [review criteria](review_criteria) to better understand how your submission will be reviewed.
```


## Checking that your paper compiles

JOSS uses Pandoc to compile papers from their Markdown form into a PDF. There are a few different ways you can test that your paper is going to compile properly for JOSS:

### GitHub Action

If you're using GitHub for your repository, you can use the [Open Journals GitHub Action](https://github.com/marketplace/actions/open-journals-pdf-generator) to automatically compile your paper each time you update your repository.

The PDF is available via the Actions tab in your project and click on the latest workflow run. The zip archive file (including the `paper.pdf`) is listed in the run's Artifacts section.

### Docker

If you have Docker installed on your local machine, you can use the same Docker Image to compile a draft of your paper locally. In the example below, the `paper.md` file is in the `paper` directory. Upon successful execution of the command, the `paper.pdf` file will be created in the same location as the `paper.md` file:

```text
docker run --rm \
    --volume $PWD/paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara
```

### Locally

The materials for the `inara` container image above are themselves open source and available in [its own repository](https://github.com/openjournals/inara). You can clone that repository and run the `inara` script locally with `make` after installing the necessary dependencies, which can be inferred from the [`Dockerfile`](https://github.com/openjournals/inara/blob/main/Dockerfile).

## Behind the scenes

Readers may wonder about the reasons behind some of the choices made for paper writing. Most often, the decisions were driven by radical pragmatism. For example, Markdown is not only nearly ubiquitous in the realms of software, but it can also be converted into many different output formats. The archiving standard for scientific articles is JATS, and the most popular publishing format is PDF. Open Journals has built its pipeline based on [pandoc](https://pandoc.org), a universal document converter that can produce both of these publishing formats as well as many more.

A common method for PDF generation is to go via LaTeX. However, support for tagging -- a requirement for accessible PDFs -- is not readily available for LaTeX. The current method used ConTeXt, to produce tagged PDF/A-3.

---

# Example Paper

This example `paper.md` is adapted from _Gala: A Python package for galactic dynamics_ by Adrian M. Price-Whelan [http://doi.org/10.21105/joss.00388](http://doi.org/10.21105/joss.00388).


```markdown
---
title: 'Gala: A Python package for galactic dynamics'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Adrian M. Price-Whelan
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Author with no affiliation
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - given-names: Ludwig
    dropping-particle: van
    surname: Beethoven
    affiliation: 3
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University, United States
   index: 1
   ror: 00hx57361
 - name: Institution Name, Country
   index: 2
 - name: Independent Researcher, Country
   index: 3
date: 13 August 2017
bibliography: paper.bib

---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# State of the field                                                                                                                  

Several tools exist for galactic dynamics computations:                                                     
`galpy` [@Bovy:2015] is a Python package with similar goals,
providing orbit integration and potential classes for galactic dynamics.                                                              
`NEMO` [@Teuben:1995] is a well-established, comprehensive stellar dynamics                                                           
toolbox written primarily in C, offering extensive functionality but with a                                                           
steeper learning curve and less integration with modern Python workflows.                                                             
Other tools like `GalPot` provide specific Milky Way potential models but lack                                                        
the broader dynamical analysis capabilities.                                                                                          
                                                                                                                                        
`Gala` was built rather than contributing to existing projects for several                                                            
reasons. First, `Gala` was designed from the ground up to integrate seamlessly                                                        
with the Astropy ecosystem, using `astropy.units` and `astropy.coordinates`                                                           
as core dependencies rather than optional features. This tight integration                                                            
enables natural workflows for astronomers already using Astropy. Second,                                                              
`Gala`'s object-oriented API with consistent interfaces across subpackages                                                            
(potentials, integrators, dynamics) provides a more modular and extensible                                                            
design than alternatives available at the time. Third, `Gala` fills a specific                                                        
niche between simple demonstration codes and full N-body simulation packages                                                          
like `Gadget` [@Springel:2005] – it focuses on the common tasks in galactic                                                             
dynamics research (orbit integration, potential evaluation, coordinate                                                                
transformations) while maintaining both performance through C implementations                                                         
and usability through its Python interface.  

# Software design

`Gala`'s design philosophy is based on three core principles: (1) to provide a
user-friendly, modular, object-oriented API, (2) to use community tools and
standards (e.g., Astropy for coordinates and units handling), and (3) to use
low-level code (C/C++/Cython) for performance while keeping the user interface
in Python. Within each of the main subpackages in `gala` (`gala.potential`,
`gala.dynamics`, `gala.integrate`, etc.), we try to maintain a consistent API
for classes and functions. For example, all potential classes share a common
base class and implement methods for computing the potential, forces, density,
and other derived quantities at given positions. This also works for
compositions of potentials (i.e., multi-component potential models), which
share the potential base class but also act as a dictionary-like container for
different potential components. As another example, all integrators implement a
common interface for numerically integrating orbits. The integrators and core
potential functions are all implemented in C without support for units, but the
Python layer handles unit conversions and prepares data to dispatch to the C
layer appropriately.Within the coordinates subpackage, we extend Astropy's
coordinate classes to add more specialized coordinate frames and
transformations that are relevant for Galactic dynamics and Milky Way research.

# Research impact statement

`Gala` has demonstrated significant research impact and grown both its user base
and contributor community since its initial release. The package has evolved
through contributions from over 18 developers beyond the original core developer
(@adrn), with community members adding new features, reporting bugs, and
suggesting new features.

While `Gala` started as a tool primarily to support the core developer's
research, it has expanded organically to support a range of applications across
domains in astrophysics related to Milky Way and galactic dynamics. The package
has been used in over 400 publications (according to Google Scholar) spanning
topics in galactic dynamics such as modeling stellar streams [@Pearson:2017],
Milky Way mass modeling, and interpreting kinematic and stellar population
trends in the Galaxy. `Gala` is integrated within the Astropy ecosystem as an
affiliated package and has built functionality that extends the widely-used
`astropy.units` and `astropy.coordinates` subpackages. `Gala`'s impact extends
beyond citations in research: Because of its focus on usability and user
interface design, `Gala` has also been incorporated into graduate-level galactic
dynamics curricula at multiple institutions.

`Gala` has been downloaded over 100,000 times from PyPI and conda-forge yearly
(or ~2,000 downloads per week) over the past few years, demonstrating a broad
and active user community. Users span career stages from graduate students to
faculty and other established researchers and represent institutions around the
world. This broad adoption and active participation validate `Gala`'s role as
core community infrastructure for galactic dynamics research.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

```

Example `paper.bib` file:

```
@article{Pearson:2017,
  	url = {http://adsabs.harvard.edu/abs/2017arXiv170304627P},
  	Archiveprefix = {arXiv},
  	Author = {{Pearson}, S. and {Price-Whelan}, A.~M. and {Johnston}, K.~V.},
  	Eprint = {1703.04627},
  	Journal = {ArXiv e-prints},
  	Keywords = {Astrophysics - Astrophysics of Galaxies},
  	Month = mar,
  	Title = {{Gaps in Globular Cluster Streams: Pal 5 and the Galactic Bar}},
  	Year = 2017
}

@book{Binney:2008,
  	url = {http://adsabs.harvard.edu/abs/2008gady.book.....B},
  	Author = {{Binney}, J. and {Tremaine}, S.},
  	Booktitle = {Galactic Dynamics: Second Edition, by James Binney and Scott Tremaine.~ISBN 978-0-691-13026-2 (HB).~Published by Princeton University Press, Princeton, NJ USA, 2008.},
  	Publisher = {Princeton University Press},
  	Title = {{Galactic Dynamics: Second Edition}},
  	Year = 2008
}

@article{gaia,
    author = {{Gaia Collaboration}},
    title = "{The Gaia mission}",
    journal = {Astronomy and Astrophysics},
    archivePrefix = "arXiv",
    eprint = {1609.04153},
    primaryClass = "astro-ph.IM",
    keywords = {space vehicles: instruments, Galaxy: structure, astrometry, parallaxes, proper motions, telescopes},
    year = 2016,
    month = nov,
    volume = 595,
    doi = {10.1051/0004-6361/201629272},
    url = {http://adsabs.harvard.edu/abs/2016A%26A...595A...1G},
}

@article{astropy,
    author = {{Astropy Collaboration}},
    title = "{Astropy: A community Python package for astronomy}",
    journal = {Astronomy and Astrophysics},
    archivePrefix = "arXiv",
    eprint = {1307.6212},
    primaryClass = "astro-ph.IM",
    keywords = {methods: data analysis, methods: miscellaneous, virtual observatory tools},
    year = 2013,
    month = oct,
    volume = 558,
    doi = {10.1051/0004-6361/201322068},
    url = {http://adsabs.harvard.edu/abs/2013A%26A...558A..33A}
}

@article{Hunt:2025,
  author = {{Hunt}, Jason A.~S. and {Vasiliev}, Eugene},
  title = {Milky Way dynamics in light of Gaia},
  journal = {New Astronomy Reviews},
  year = 2025,
  volume = 98,
  doi = {10.1016/j.newar.2024.101721},
  url = {https://www.sciencedirect.com/science/article/pii/S1387647324000289}
}

@misc{fidgit,
  author = {A. M. Smith and K. Thaney and M. Hahnel},
  title = {Fidgit: An ungodly union of GitHub and Figshare},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {https://github.com/arfon/fidgit}
}
```

Note that the paper begins by a metadata section (the enclosing --- lines are mandatory) and ends with a References heading, and the references are built automatically from the content in the `.bib` file. You should enter in-text citations in the paper body following correct [Markdown citation syntax](https://pandoc.org/MANUAL.html#extension-citations).  Also note that the references include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
