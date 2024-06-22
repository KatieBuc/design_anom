<h1>Design and Anomaly Detection in River Networks</h1>

<h2>Overview</h2>
<p>This repository contains the code and data for the paper "Bayesian Design for Sampling Anomalous Spatio-Temporal Data" <a href="#footnote1" id="ref1">[1]</a>. The primary objective of this project is to develop a robust Bayesian optimal experimental design (BOED) framework with anomaly detection methods to ensure high-quality data collection. The framework involves anomaly generation, detection, and error scoring within the search for optimal designs. See below for implementation for two simulated case studies: a spatial dataset and a spatio-temporal river network dataset.</p>

<h2>Repository Structure</h2>
<pre>
.
├── river_spatiotemporal_anomaly_design.R
├── spatial_anomaly_design.R
└── README.md
</pre>

<h2>Getting Started</h2>

<h3>Prerequisites</h3>
<p>Ensure you have the following software installed:</p>
<ul>
    <li>R (version 3.0 or higher)</li>
    <li>RStudio (optional but recommended)</li>
    <li>Required R packages: <code>oddstream</code>, <code>dplyr</code>, <code>dplyr</code>, <code>SSN</code>, <code>FNN</code></li>
</ul>

<h3>Installation</h3>
<ol>
    <li>Clone this repository to your local machine:
        <pre><code>git clone https://github.com/KatieBuc/design_anom.git
cd design_anom
        </code></pre>
    </li>
    <li>Install the required R packages. Open R or RStudio and run:
        <pre><code>install.packages(c("oddstream", "dplyr", "dplyr", "SSN", "FNN"))</code></pre>
    </li>
</ol>

<h3>Usage</h3>
<p>The scripts are designed to be run on a High-Performance Computing (HPC) cluster. Below are the steps to submit the scripts to the HPC.</p>

<strong>River network example</strong>: 
        <pre><code>R -e "seed <- $seed; prop=$prop; lambda=$lambda; scenario=$scenario; source('./river_spatiotemporal_anomaly_design.R');</code></pre>


<h2>Contributing</h2>
<p>If you would like to contribute to this project, please fork the repository and create a pull request with your changes. Contributions are welcome and appreciated!</p>

<h2>License</h2>
<p>This project is licensed under the MIT License. See the <code>LICENSE</code> file for more details.</p>


<p>Thank you for using this repository! We hope you find it useful for your research and projects.</p>

<p id="footnote1"><sup>1</sup> Buchhorn, Katie, et al. <i>"Bayesian Design for Sampling Anomalous Spatio-Temporal Data."</i> arXiv preprint arXiv:2403.10791 (2024). <a href="#ref1">↩</a></p>

[^1]: Buchhorn, Katie, et al. _"Bayesian Design for Sampling Anomalous Spatio-Temporal Data."_ arXiv preprint arXiv:2403.10791 (2024).
