# Preface

This textbook is about **PlaSiC**, a numerical general circulation model I developed for teaching. PlaSiC is a model of moderate complexity: it consists of a complete codebase of about 20,000 lines of C, accompanied by this textbook, which explains in detail the principles behind the code.

This naturally raises a question: why build a model like this from scratch at all?

Developing PlaSiC took a great deal of work. People often tell me this is a thankless undertaking. Research in China is fiercely competitive, and everyone is racing to publish. Why pour so much time into something that is meant mainly for teaching (and so offers little novelty) and is hard to turn into strong publications (and so brings little immediate return)?

There are several reasons.

## 1. A regret I could never shake

One regret has stayed with me for years: I never learned how to build a real GCM from scratch.

I did my undergraduate degree at the School of Atmospheric Sciences, Nanjing University, where great emphasis is placed on solid foundations in mathematics and physics, as well as on numerical modeling. I was fortunate to receive a rigorous, systematic education in atmospheric science there — yet it was precisely this training that made me increasingly aware of the shortcomings in how numerical models are taught. Let me start with my own story.

My first encounter with a model came in the second semester of my sophomore year, in the first programming course I ever took: *C Programming* (《C语言程序设计》). The final assignment was a group project — write a program in C and present it to the class. Most groups built something like a Snake game or a student information management system. My group — Lian Nayuan (连娜媛), Liu Qi (刘淇), and I — chose instead to build a very simple model: a "Daisyworld" Earth system simulation. I remain grateful to them both. Life brought us together and then carried us in different directions, and we have long since lost touch.

Daisyworld is an imaginary planet covered with black and white daisies. The black daisies absorb sunlight and warm the planet, while the white ones reflect it and cool the planet, so the shifting balance between them keeps the global temperature remarkably stable. The model was devised as a simple demonstration of the **Gaia hypothesis**, proposed by James Lovelock in the 1960s. The hypothesis holds that the Earth behaves in some ways like a living organism: life and its environment interact with and regulate each other, together forming a stable, self-regulating system. And the more complex the system — that is, the greater the planet's biodiversity — the better it can withstand external disturbances.

<div class="figure-row figure-row--two">

  <figure>

    <img src="../assets/images/daisy_1.png" alt="Title page and design notes of the Daisyworld Earth system simulation report">

  </figure>

  <figure>

    <img src="../assets/images/daisy_2.png" alt="Data structures and screenshots of the Daisyworld game interface">

  </figure>

  <p class="figure-row__caption">Fig. 0.1.1. <strong>Our Daisyworld model.</strong> Left: the design report for our final project in the C programming course. Right: screenshots of the program's interface — the start screen, the dialog for setting the initial daisy coverage, and the year-by-year simulation display.</p>

</div>

I would especially like to thank Prof. Wang Shuyu (王淑瑜) of the School of Atmospheric Sciences. She often sat down with us to discuss how to implement the simulation, effectively serving as our scientific advisor. Without her help, we could never have finished the project.

Naive as it looks to me now, that project planted the first seed of an idea in my mind: simulating the Earth.

Later, in the first semester of my junior year, I began to study numerical models more systematically in *Numerical Weather Prediction* (《数值天气预报》), a course taught by Prof. Sun Xuguang (孙旭光). I worked hard all semester, but to be honest, I only half understood the material — and by the end of the following summer, I had forgotten almost all of it.

Still unwilling to give up, I took another course in the autumn of 2022, the first semester of my senior year: *Atmospheric Numerical Simulation Experiments* (《大气数值模拟实验》), taught by Yang Ben (杨犇). The course focused mainly on running numerical experiments with WRF (the Weather Research and Forecasting model) rather than on the principles behind numerical models.

In one lecture, the instructor mentioned an idealized WRF experiment in which entire mountain ranges had been flattened. I found this fascinating, but at the time I could not imagine how such an experiment was actually carried out — although it seems straightforward enough to me now. Sadly, even after finishing the course, I still understood very little about how numerical models actually work.

In 2023, I went straight from my undergraduate studies into a PhD program in the Department of Earth System Science at Tsinghua University, where I took yet another course: *Earth System Numerical Modeling* (《地球系统数值模式》). It was taught in relay fashion, with each lecture given by a different professor on their own area of expertise. On paper, this "all-star lineup" sounded impressive. In practice, however, I found the course disappointing — I would go so far as to call it a failure.
Each lecture gave us another fragment of knowledge, but nothing ever tied the fragments together into a coherent framework; for me, at least, they never added up to a systematic understanding. At times I even suspected that many of the instructors had been away from hands-on model development for so long that they themselves no longer fully understood its practical details.

In short, I had come a long way and still did not truly understand numerical models. To me, they remained a black box: I had learned how to run one, but never what went on inside.

I doubt I am not alone in this. In China, and indeed around the world, there must be countless students like me: they take the courses, come away with only a hazy understanding of models, and end up intimidated by them. After carefully surveying the available teaching resources, I concluded that there is a real gap: good teaching models — ones that come with complete source code and a complete tutorial, and that are small enough to understand yet include every essential component — are extremely rare.

Jumping straight into a mature system such as CESM2 (the Community Earth System Model, version 2) is unrealistic. With millions of lines of code, it is far too large for any student to study and fully understand in a few months. A mature, complex model must first be distilled into something students can actually digest.

So I decided to build one myself.

## 2. Helping machine-learning researchers enter AI for weather and climate

My second motivation comes from the rapid progress of machine learning in recent years, which has opened up a new research field: AI for weather and climate (AIWC).

Today, a substantial share of AIWC researchers have no formal training in atmospheric science; many come from computer science and related fields instead. Their programming skills and command of AI are often excellent, but they may lack sufficient domain knowledge — in particular, many are unfamiliar with the workflow and conceptual framework of conventional numerical weather prediction. This gap is one reason why some AIWC papers look baffling to meteorologists.

I believe one of the fastest ways to close that gap is to build an NWP model yourself, from the ground up.

I therefore hope that PlaSiC can serve as a teaching tool for the wider machine-learning community, helping more researchers enter AIWC with a firmer grasp of the atmospheric science that underpins it. In this way, I hope it can contribute, however modestly, to the growth of the AIWC community as a whole.

## 3. A bitter lesson from my own research

My third motivation comes from my own research. Since starting my PhD in August 2023, I have spent three years working on AIWC.

I once hoped to build a fully automated, AI-driven operational system covering the entire chain, from data assimilation to weather, subseasonal-to-seasonal, and climate prediction (see [About the Author](author.md) for details). As my research progressed, however, I gradually came to a bitter realization: future AIWC systems may still need a numerical model at their foundation — or, more precisely, as their dynamical core.

I believe PlaSiC is a useful starting point for exploring this possibility, and in the future I may integrate it more closely with AI. That, however, lies beyond the scope of this textbook, so I will not dwell on my research plans here.

## 4. A pragmatic reason

My fourth reason is, frankly, a pragmatic one: I need a substantial open-source project that shows what I can do.

I have interned at Caiyun Weather and at Huawei's 2012 Laboratories, working on spatiotemporal prediction. Those experiences taught me a simple truth: companies ultimately care about creating value, so they need the people they hire to deliver useful results within a realistic timeframe. In other words, what matters is having the skills that companies actually need.

Even a paper in *Nature* may not carry much weight if it is a conventional climate-analysis study along the lines of "under global warming, X changes, which in turn has adverse impact Y on human society." Such work may be valuable as academic research, but it does not necessarily demonstrate the capabilities an employer is looking for.

At the same time, as AI agents spread rapidly and grow ever more capable, the traditional way of evaluating programmers — having them write code by hand — is becoming obsolete. We used to be tested on algorithmic programming problems, but many of those problems are now trivial for modern AI systems.

I believe that in the future, experience in designing and maintaining large software systems, managing complex projects, and making effective use of large-scale computing resources will become increasingly valuable. For these reasons, I hope PlaSiC will also be a meaningful addition to my CV.

## 5. Thanks to Prof. Jiang Yanyan

I owe special thanks to Prof. Jiang Yanyan (蒋炎岩) of the School of Software at Nanjing University. Prof. Jiang is a passionate teacher who has posted many excellent lecture videos on Bilibili, including those for his operating systems course.

For that course, he did things the hard way: he wrote more than 100,000 lines of code by hand to build a framework in which undergraduates implement their own compiler, operating system, and processor — and then run their own application ecosystem on top of them.

Prof. Jiang has said that, by the standards of conventional academic research, this is about the least cost-effective thing one could possibly do — after all, teaching is often undervalued in Chinese universities. Yet to him, it is the most worthwhile investment of all.

His attitude inspired me deeply and became an important motivation for developing PlaSiC.

![Prof. Jiang Yanyan speaking about building the operating systems course from scratch](assets/images/jiangyanyan.png)

Fig. 0.1.2. **"The least cost-effective thing one can do" — and yet the most cost-effective.** Prof. Jiang Yanyan (School of Software, Nanjing University), speaking in the talk series *All About Doing a PhD* (《读博那些事儿》). For his operating systems course, he wrote more than 100,000 lines of code to build a framework in which undergraduates implement their own compiler, operating system, and processor, and then run their own application ecosystem on top of them.

## Closing words

These are the reasons I built PlaSiC. The model itself is introduced in detail in the [PlaSiC Online Tutorial](index.md), so I will leave the technical details there.

Finally, I would like to thank the teachers who helped me along the way — and to thank myself, too, for piecing this project together over the past nine months. 
I hope PlaSiC will prove useful to a wider community.

Thank you for reading.
