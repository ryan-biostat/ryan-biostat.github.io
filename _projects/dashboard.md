---
layout: page
title: Lab Sample Dashboard
description: a Python + JavaScript dashboard for tracking whole genome samples as they process from lab to analysis
img: assets/img/charts.jpg
importance: 2
category: Biostatistics
---

#### **About**

To improve lab communication, we required a process to track samples processed for Long-Read Sequencing (LRS) (see [LRS Pipeline](Ryan-biostat.github.io/projects/ONT_pipeline/)).  My thought was that an easily accessible, browser-based dashboard so anyone in the lab could quickly check sample status without digging through files or logs. The only problem was that I hadn't had much familiarity with tools that would create dashboards.  I also didn't have the expertise to design one from scratch or host it locally. 

This presented a unique opportunity to learn a new language through AI-assisted development. With heavy support from Google Gemini, I built a Vite + React dashboard that pulls data via APIs and Python scripts to retrieve, clean, and format the underlying data. The result is a friendly, easily-navigable dashboard that communicates LRS processing that give me hands-on experience with webpage development and LLM-assisted coding.

## **Purpose**

Our lab currently manages nearly a hundred samples at various stages of the Long-Read Sequencing workflow, including completed samples, those in progress, and others still waiting to be processed. Clinicians and collaborating physicians often request updates, and answering these questions requires navigating our Linux-based computing cluster, which is not accessible to everyone. This dashboard was created to give any member of the lab a clear and straightforward view of the information available for each sample. It addresses our immediate communication needs while providing a flexible foundation for future growth, allowing new metrics and visualizations to be added as the lab’s questions and workflows evolve.



## **Methods**

1. #### **Figure Out What to Build**

   Before writing any code, I spent time exploring how best to create an interactive and shareable dashboard for the lab. I initially looked at business intelligence tools (PowereBI, Tableau) but they were slow and costly. Eventually, this search led me to using Vite and React that was advertised to be fast and flexible. Since I was new to JavaScript and React, I knew I'd have to rely heavily on large language models (LLMs). These tools offer help designing components, troubleshooting errors, and translating ideas into working code.

2. #### **Select an LLM**

   I compared several models focusing  coding, responsiveness, and debugging efficiency. Based on benchmarks and hands-on testing, I chose Gemini Pro 2.5, which consistently outperformed competitors at the time. Another major advantage was its Canvas mode that lets users edit code directly, highlight specific sections, and ask targeted questions within the same interface. This made it easy to iterate quickly and fix errors without jumping between multiple tools.

3. #### **Generate a Dashboard and Customize**

   I connected the project’s GitHub repository directly to Gemini, which made it easy for the model to read the code, propose edits, and build new features. Reaching the design I envisioned took a lot of back-and-forth until everything behaved the way I wanted. Through that iterative process, I ended up learning quite a bit of JavaScript and gaining a clearer understanding of how React components work together. The final structure of the dashboard came from this cycle of generating code, reviewing the results, and gradually shaping it into something polished and functional.

4. #### **Gather Data from Local Sources**

   All the information displayed on the dashboard came directly from our lab's Excel sheets. To make this usable in a web dashboard, I wrote Python scripts to parse multiple spreadsheets, clean the fields, and aggregate everything into a single structured document that the dashboard could easily read. Because new samples are constantly added and statuses change over time, this workflow needed to be reliable and repeatable. The scripts were designed to run the same process over and over, pulling fresh data each time without breaking, so the dashboard would always reflect the most current information available.

5. #### **Publish and Maintain**

   Once the dashboard was functional, the final step was making it accessible to the rest of the lab. I deployed it through our computing cluster so it could be viewed on our local network without exposing anything externally. Since the dashboard displays sensitive sample information, I worked closely with our security team to ensure that all access was restricted and that no data could be reached from outside our internal environment. I also set up the Python retrieval script to run automatically every morning, pulling in the latest Excel updates and refreshing the dashboard’s data. This automated process has been running reliably for nearly six months, keeping the dashboard current without requiring any manual intervention.

## **Snapshot**

Here is a snapshot of our project dashboard. Each sample populates a single row with icons communicating their status in the processing pipeline.

{% include figure.liquid 
    path="assets/img/dashboard_snap.png"
    title="dashboard Image"
    class="img-fluid rounded z-depth-1"
%}



