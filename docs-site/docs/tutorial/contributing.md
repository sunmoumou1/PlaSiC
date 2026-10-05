# How to Vibe-Code Your Own GCM — Lessons from Experience

This document is for people who want to vibe-code their own GCM, or who want to make major additions or improvements to PlaSiC. If all you want to do is use PlaSiC to run experiments, you can safely skip this document.

With coding agents exploding in capability and popularity since 2026, the traditional scientific workflow is changing at remarkable speed. I think 2026 will be remembered as a watershed year: from this point on, research will rely heavily on AI agents, especially for coding. In the future, master's and PhD students—the people who do much of the day-to-day work of research—may hardly write code by hand at all. Quite possibly, not a single line.

That means we need to learn how to work with AI agents as collaborators, and how to use them sensibly in research and creative work. From 2026 onward, work you develop with an AI agent is still your own work and your own creation. Ten years ago, nobody would have argued that code written with an IDE's autocomplete somehow belonged to the IDE. AI agents are far more capable, of course, but the basic principle is not so different.

I will not go into every detail here. Instead, I will focus on a few basic principles—and these principles matter a great deal. In fact, I think they matter more than any step-by-step recipe for development. The exact tools and workflows will change quickly; the principles probably will not.

## Principle 1

!!! warning "⚠️ GCMs are extremely complex! ×3"

!!! info "🧠 So you must respect that complexity in your head! ×3"

!!! danger "🙋 So you must keep a human in the loop! ×3"

A good analogy is an operating systems course in a computer science curriculum. Operating systems is notoriously difficult, but learning it pays off enormously later on. **It marks the transition from building small, self-contained programs to thinking seriously about how large, complex software systems actually work.**

Processes, threads, file systems, virtual memory, disk I/O, and all the other concepts you encounter are like the individual parts of a car. What an operating systems course teaches you is how those parts fit together into an actual machine—one that not only runs, but runs reliably.

The same is true in atmospheric science. Building a GCM yourself—or, more broadly, learning how a GCM works and how it is implemented—is as important to an atmospheric scientist as studying operating systems is to a computer scientist.

Only after building a GCM yourself do you begin to see just how extraordinarily complex the atmosphere and climate system really are. A student may read a handful of papers and come to understand one particular physical process—for example, the formation and propagation of midlatitude Rossby waves. After studying it long enough, they may even feel that the process is fairly well understood and, in some sense, under control.

But the atmosphere contains an enormous number of such processes. No single person has enough time or mental bandwidth to understand every one of them in full detail, much less to keep track of all the ways they interact and couple with one another.

And this is precisely why building a large, complex scientific program is not something you can accomplish by throwing a few prompts at Claude Code or Codex. At least, not in 2026—and personally, I do not expect coding agents to reach that point within the next few years either.

So you need to keep a human in the loop.

<p align="center">
  Human in the loop!
</p>

<p align="center">
  Human in the loop!
</p>

<p align="center">
  Human in the loop!
</p>

Before development begins, you need a clear mental model of the software itself: what components it contains, what role each component plays in the system as a whole, in what order those components are called, how they interface with one another, and what variables or data structures are passed between them.

Thinking these questions through in advance is crucial.

For example, before you start development, you should already have a framework something like this in your head:

![Experiment design and prescribed CO₂ forcing](./assets/images/plasic_atmospheric_parameterizations.png)

Fig. 0.3.1. **Call order of PlaSiC.** Call order of the dynamical core and the physical parameterizations in one PlaSiC time step, with the variables exchanged between them.



## Principle 2

A coding agent is an assistant, not Aladdin’s lamp. How far you can go with an agent still depends, to a large extent, on how far you could get on your own.

If you could eventually complete, say, 60% of a project by yourself—given enough time and effort—then a coding agent can give you wings. But if you would have absolutely no idea where to begin without one, then chances are you will not get very far with one either.

So before anything else, you need to understand how a GCM actually works.

Here I strongly recommend four resources. In fact, it was only after working through these four carefully that I personally felt ready to take on the design and development of a GCM:


![Experiment design and prescribed CO₂ forcing](./assets/images/resources.png)

Fig. 0.3.2. **Recommended resources.**

If you want to build your own GCM from scratch, I strongly recommend starting with the four resources above. At the very least, work through one of them in full—or use other resources that you find more suitable.

## Principle 3

**Being the engineer of a project is not only writing code; I can say plainly that writing code takes only a small part of a project's development and maintenance time. An engineer's main responsibility is to design, develop, maintain, apply and promote a project. An engineer is responsible for the entire lifecycle of the project — responsible! None of these steps is simple.**

So before you try to develop your own GCM, be fully mentally prepared: this is a hard, protracted war. (In fact, I developed PlaSiC from February 2026 until the end of September 2026.) And make sure you can coordinate a large amount of sustained time investment — a large amount of sustained time investment, a large amount of sustained time investment.

## Principle 4

!!! info "⚙️ To use a tool well, you must understand how it works! ×3"

If we want to use coding agents to help us develop a GCM, we must at least understand roughly how coding agents are implemented (no need for fine detail).

The essence of a coding agent is:

!!! info "Agent = LLM + Context + Tools"

Understanding this essence will greatly improve your efficiency in using coding agent tools! (Very important.)

This formula comes from *Understanding AI Agents: Design Principles and Engineering Practice* by Bojie Li ([github.com/bojieli/ai-agent-book](https://github.com/bojieli/ai-agent-book)). I strongly recommend that project, or other lectures and videos; in short, it is best if you understand how LLMs work and how agents work.

Here is one concrete example of why this matters: context is crucial.

Suppose you ask a coding agent to implement the Kuo convection parameterization. You could place a `kuo_convection.md` file in the repository containing clear notes you made while reading the classic papers. You could also create a `reference_code/` directory with implementations from SpeedyWeather.jl or other models.
That way, the coding agent has immediate access to the relevant theory, your own understanding of it, and working reference implementations it can inspect while developing the new module.

## Principle 5

!!! tip "🎯 In early development, only results matter! ×3"

Understanding this can make your development process dramatically more efficient. For example:

1. When implementing a new feature, do not obsess over the exact code at the beginning. What matters first is that you understand, at a high level, what kind of approach could solve the problem. Understanding the underlying principle and writing the actual implementation are two different things. Instead, spend your effort designing good test cases and clear acceptance criteria. Then let the agent—or a whole swarm of agents—work aggressively toward the goal until they believe the task is complete and every test or requirement you defined has been satisfied.

2. You can go a step further by bringing in other agents as reviewers. If agent *a*, powered by LLM A, writes the code, then have agent *b*, powered by LLM B, review it independently and verify that the implementation is actually complete and correct. You can even assign the same task to several agents independently and keep only the best implementation in the Git repository. I used exactly this strategy when testing PlaSiC's weather-forecasting capability. I had three agents, each powered by a different LLM, independently build the same end-to-end workflow: automatically download ERA5 data, preprocess it, perform data assimilation, run the forecast, evaluate the results, and generate plots. In the end, I kept the implementation that produced the smallest forecast error.

## Principle 6

**Review code selectively and in layers.**

I strongly advise against reviewing a codebase by reading every file line by line from beginning to end. But the opposite extreme—reading no code at all and leaving everything to the agent—is even worse. That is how a codebase turns into a big ball of mud surprisingly quickly.

Review with clear priorities. Start with the core logic. Simple utilities, boilerplate, and non-critical code can often be skimmed, or sometimes skipped entirely.

More importantly, review strategically. Do not begin by diving straight into the raw source code. First, ask the agent to produce a `pseudocode.md` for the part you care about: a natural-language explanation of the implementation, including the pseudocode, algorithmic flow, key data transformations, and an overall walkthrough.

This gives you two major advantages. First, you can build a high-level mental model of the implementation very quickly, which makes the detailed code review much faster afterward. Second, conceptual or architectural flaws often become obvious at this stage, before you spend time digging through implementation details that may need to be rewritten anyway.

## Principle 7

Tokens are not free, so save them where you can. There is no reason to throw your most capable—and most expensive—LLM at every task. Learn to match the model to the difficulty of the problem.

My day-to-day development stack looks roughly like this. For difficult planning and architectural work, I use GPT-6 Astra. For complex code modifications, debugging, and testing, I use GPT-6.1 Sol, usually with the reasoning level set to Ultra or even Extra High. For routine development tasks—which make up the vast majority of my workload—I use DeepSeek V4.1 Flash. And for very simple questions, I use the cheapest option available, GPT-6 Luna with Extra High reasoning.

**In other words: do not use a sledgehammer to crack every nut.**

Of course, knowing which model is appropriate for which task is not something you learn overnight. It takes a lot of hands-on use before you develop a reliable intuition for when a stronger model is genuinely worth the extra cost—and when it is simply overkill.

## Principle 8

**Review the entire codebase periodically.** This is important—important enough that it can determine whether a project actually makes it to completion or gets abandoned halfway through.

For example, set aside time once a month to review the repository as a whole. During these reviews, focus less on individual lines of code and more on strengthening your mental model of the project: its architecture, major modules, dependencies, data flow, and control flow.

The ultimate goal is to reach the point where you can mentally unfold a map of the entire codebase.

<p align="center">
  The whole repository! In your head!
</p>

Regularly clean up duplication as well. One common weakness of coding agents is their tendency to generate too much defensive code, along with even entire files that are unnecessary. Remove this kind of redundancy as soon as you notice it.

Always remember: every extra line of code increases the long-term maintenance burden. A healthy codebase is not one that only grows—it is one that is regularly pruned.

## Principle 9

!!! tip "👣 Get your hands dirty — a journey of a thousand miles begins with a single step! ×3"

!!! note "🧩 Features are not completed in one step; they are added bit by bit! ×3"

!!! info "🚀 Get going, don't fear failure — just charge ahead! ×3"
