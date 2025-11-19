---
layout: page
permalink: /teaching/
title: teaching
description: Materials for courses taught.
nav: true
nav_order: 3
---

To achieve the look where every course is its own distinct "container" with a header (title/subtitle) followed by a list of materials, the best approach in al-folio (which uses Bootstrap) is to use the Bootstrap Card component.

Here is the code to produce a "Teaching" page with a boxed layout for each course.

_pages/teaching.md
Copy and paste this entire block into your file.

HTML

---
layout: page
title: teaching
permalink: /teaching/
description: materials for courses you taught.
nav: true
nav_order: 5
---

<h3>Carnegie Mellon University</h3>

<div class="card mt-3 mb-3">
  <div class="card-header">
    <div class="row align-items-center">
      <div class="col-md-6">
        <strong>Probabilistic Graphical Models</strong>
      </div>
      <div class="col-md-4 text-muted">
        Spring 2020: Guest Lecturer
      </div>
      <div class="col-md-2 text-right text-muted">
        10-708
      </div>
    </div>
  </div>
  
  <div class="card-body p-0">
    <table class="table table-sm table-borderless mb-0">
      <tbody>
        <tr>
          <td class="pl-4">Lecture 19: Reinforcement Learning (Part 1)</td>
          <td class="text-right pr-4">
            <a href="#" class="badge bg-primary">slides</a>
            <a href="#" class="badge bg-secondary">video</a>
            <a href="#" class="badge bg-success">notes</a>
          </td>
        </tr>
        <tr>
          <td class="pl-4">Lecture 20: Reinforcement Learning (Part 2)</td>
          <td class="text-right pr-4">
            <a href="#" class="badge bg-primary">slides</a>
            <a href="#" class="badge bg-secondary">video</a>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</div>

<div class="card mt-3 mb-3">
  <div class="card-header">
    <div class="row align-items-center">
      <div class="col-md-6">
        <strong>Intro to Machine Learning</strong>
      </div>
      <div class="col-md-4 text-muted">
        Spring 2019: Head TA
      </div>
      <div class="col-md-2 text-right text-muted">
        10-701
      </div>
    </div>
  </div>

  <div class="card-body p-0">
    <table class="table table-sm table-borderless mb-0">
      <tbody>
        <tr>
          <td class="pl-4">Recitation: Neural Networks</td>
          <td class="text-right pr-4">
            <a href="#" class="badge bg-dark">code</a>
            <a href="#" class="badge bg-primary">slides</a>
          </td>
        </tr>
        <tr>
          <td class="pl-4">Homework 3: SVMs</td>
          <td class="text-right pr-4">
            <a href="#" class="badge bg-info">solution</a>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</div>