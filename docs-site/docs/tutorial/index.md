<p class="tutorial-logo">
  <img src="assets/images/logo.png" alt="PlaSiC logo">
</p>

# PlaSiC Online Tutorial

!!! tip "Start here"
    New to the project? Read the [Preface](preface.md) first: it tells the story of why PlaSiC was hand-built and what this tutorial is trying to fix.

## Who this tutorial is for

The tutorial addresses two kinds of readers:

- The first group consists of senior undergraduate students majoring in atmospheric science, as well as graduate students working on the development of numerical weather and climate models. In many university curricula, students have limited access to well-designed courses on numerical modeling. Although they may learn some basic numerical methods, they often find the material abstract and confusing, forget much of it soon after completing the course, and still lack a clear understanding of how climate models are actually implemented from an engineering perspective. At least, this was my own experience as a student.
For these learners, this tutorial provides an accessible entry point. It includes a complete numerical model with a manageable codebase of approximately 20,000 lines, making it feasible to study the model as a whole, together with accompanying online learning materials.

- The second group consists of machine learning researchers who plan to work on ML for climate. These researchers may want to develop data-driven AI models for applications such as weather forecasting and climate simulation. While they are often highly skilled in machine learning and computer science, they may lack sufficient domain knowledge in atmospheric science and numerical modeling.
PlaSiC is therefore a useful starting point for this audience as well. It can help ML researchers quickly become familiar with numerical model development while developing a deeper understanding of how numerical weather and climate models are implemented in practice.

![The two kinds of readers this tutorial addresses: climate-science students and machine-learning researchers entering AI for weather and climate](assets/images/user_portrait.png)

Fig. 0.2.1. **The two kinds of readers this tutorial addresses.** The first group, on the left, wants to understand how a GCM is actually built rather than only how to run one; the second, on the right, enters AI for weather and climate (AIWC) and needs the domain background that only a complete numerical model can provide.

## Key features of PlaSiC

- **A model of moderate complexity.** PlaSiC is a fully coupled climate system model written in about 20,000 lines of code. It is detailed enough to reproduce the essential behaviour of the climate system, yet compact enough for a single reader to study from beginning to end. The tutorial therefore follows the model module by module, from the spectral dynamical core to the physical parameterizations, in the exact order in which they are called within one time step.

    ![Call order of the dynamical core and the physical parameterizations within one PlaSiC time step](assets/images/plasic_atmospheric_parameterizations.png)

    Fig. 0.2.2. **Structure of a PlaSiC time step.** We will cover the all the modules in detail in the next chapters.

    ![Word cloud of the governing equations used by the PlaSiC dynamical core and physical parameterizations](assets/images/equations_word_cloud.png)

    Fig. 0.2.3. **The mathematics behind the model.** A word cloud of the governing equations covered by this tutorial. There are a total of thirteen equations and thirteen unknowns, forming a closed equation system.

- **A complete and detailed textbook.** These online materials are more than a user manual: every component is derived from its governing equations and then mapped onto the corresponding source files, so that the reader learns both the physics and its implementation. Note that this textbook will be updated as the model evolves, so it is not a static document.

- **An interactive visual interface.** The model is shipped with PlaSiC Studio, a graphical front end that lets you set the horizontal resolution, vertical levels, integration length and initial data, build the model, and start or stop a run with a click. Results are streamed while the model is running and can be inspected as 2-D fields on a rotatable globe or as 3-D volumes over a selected region, which makes it easy to build intuition before working with the raw model output.

    ![PlaSiC Studio: experiment configuration on the left and the live climate field on an interactive globe on the right](assets/images/app1.png)

    Fig. 0.2.4. **PlaSiC Studio: configuration and live 2-D visualization.** The left panel controls the build and run; the right panel shows a selected field, here wind speed, evolving on the globe.

    ![PlaSiC Studio: 3-D volume rendering of the wind field over a selected longitude-latitude-sigma region](assets/images/app2.png)

    Fig. 0.2.5. **PlaSiC Studio: 3-D visualization.** The same run can be explored as a volume rendering over a chosen region and vertical extent.

- **CMIP-style experiments out of the box.** PlaSiC can be used to carry out standard climate experiments such as the historical 1850–2014 simulation, including the spin-up and control runs required to interpret the forced response. The experiment scripts are provided, so that each experiment documented in [Chapter 6](6-experiments/6.2-historical-experiment-1850-2014.md) can be reproduced with a single command.

## How PlaSiC compares with other models

Several excellent systems occupy the space between teaching and research climate modelling. The table below summarizes how PlaSiC relates to the best-known of them. The comparison is deliberately opinionated: its purpose is to help you choose the right tool for your goal, not to rank the projects.

| Model | Original role | Model scope | Workflow | PlaSiC perspective |
| --- | --- | --- | --- | --- |
| **PUMA** | University teaching in atmospheric dynamics and numerical methods. | 3-D primitive equations with Newtonian cooling and Rayleigh friction; primarily a dry dynamical core. | Fortran with a traditional GUI. | Excellent for teaching primitive-equation dynamics, but with a lighter treatment of physical parameterizations and climate-system coupling. |
| **PlaSim** | Training in general circulation model development. | A primitive-equation AGCM with radiation, clouds, convection, land, a mixed-layer ocean, and sea ice. | Fortran with a GUI. | PlaSiC rebuilds this foundation in C, making the source more approachable to ML researchers and pairing it with detailed teaching material and a modern visualization app. |
| **Isca** | A climate-model hierarchy for graduate teaching and research. | Ranges from Held–Suarez experiments to moist GCMs with radiation, land, and a mixed-layer ocean. | Fortran core with Python tooling. | Its hierarchy is powerful for comparing model complexity. PlaSiC instead focuses on mastering one coherent, complete GCM. |
| **CliMT** | A climate-modelling framework for students and researchers. | Independent components for dynamics, radiation, convection, surface processes, ice, and snow. | Python API and Jupyter. | CliMT excels as a modular construction kit. PlaSiC is designed to be read end to end—from governing equations and numerics to physics and coupled surfaces. |
| **climlab** | Undergraduate and graduate teaching in climate physics. | Energy-balance, radiative–convective, and 1-D/2-D process models rather than a full 3-D primitive-equation GCM. | Python and Jupyter. | A mature choice for explaining climate processes with low-dimensional models; PlaSiC reveals how those ideas operate inside a complete 3-D GCM. |
| **SpeedyWeather.jl** | Modern, interactive atmospheric and climate modelling for teaching and research. | Barotropic, shallow-water, and dry/moist primitive-equation models with convection, clouds, radiation, land, slab ocean, and sea ice. | Julia, REPL, and Pluto; highly modular. | The closest modern alternative in scope. PlaSiC differentiates itself through a C implementation and a dedicated graphical application. |
| **EdGCM** | Climate-change experiments with a real GCM for students. | Built on NASA GISS Model II rather than redesigned specifically for source-level comprehension. | Desktop GUI; now legacy software. | EdGCM lowers the barrier to running a professional GCM. PlaSiC prioritizes understanding and modifying the model itself. |
| **EzGCM** | Cloud-based Earth-system modelling for education. | A web platform for running climate models and analysing their output. | Web GUI and cloud execution. | EzGCM minimizes installation and computing friction, with emphasis on simulation and visualization. PlaSiC extends the learning path into equations and source code. |

Among these, **SpeedyWeather.jl** is the closest modern peer in terms of scope, and PlaSiC is complementary to it. SpeedyWeather.jl leans towards weather prediction and places strong emphasis on simulation speed and efficiency. PlaSiC, by contrast, places greater emphasis on understandability: it is written in C, a language with a much broader global audience, and its workflows lean towards climate—it includes a three-dimensional ocean and scripts for reproducing the historical experiment. In short, readers are warmly encouraged to study SpeedyWeather.jl as supplementary material.

## Planned features and roadmap

PlaSiC is an actively developed model. The following features are not yet implemented, but are planned for future releases:

- **Aerosol-related physical processes.** 
- **External ozone forcing.** Ozone will be changed from an internally prescribed distribution to a field that can be read from external data, so that historical experiments can be driven by both ozone and carbon dioxide forcing simultaneously.
- **Particle tracking and visualization.** A particle tracking algorithm will be added, together with visualization support, so that the transport of air parcels and tracers can be followed through the simulated circulation.
- **Weather simulation experiments.** A forecast-mode experiment will be added in which the model is initialized from an observed weather state and integrated forward. The growth of forecast error as a function of lead time, from day 0 to day 10, can then be quantified.

## How to report problems

- If you find documentation errors, broken links or formula rendering problems, please open an issue on [GitHub Issues](https://github.com/sunmoumou1/plasic/issues);
- If you want to improve the documentation or tutorial content, pull requests are welcome;
- For problems related to running the model, please attach `run.log` and the relevant configuration to the issue so that it can be diagnosed.

