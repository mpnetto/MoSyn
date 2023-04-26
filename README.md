# MoSyn: Time-Varying Graph Analysis Tool

MoSyn is a MATLAB-based application designed for the analysis of time-varying graphs (TVGs) and their associated graph measures. The tool provides a modular structure with various classes and functions to handle different aspects of the analysis, such as configuration, graph features, and project management.

## Features

- Read and process time-varying graph data from input files.
- Calculate various graph features and measures, such as adjacency spectral norms (ASN) and others registered in the `MeasureRegistry`.
- Organize the data into a project structure with groups and subjects.
- Write the calculated graph features and measures to output files.
- Graphical User Interface (GUI) for easy user interaction.

## Getting Started

### Prerequisites

- MATLAB (version R2016b or later is recommended)

### Installation

1. Clone the repository to your local machine or download it as a ZIP file and extract it.
```
git clone https://github.com/yourusername/mosyn.git
```

2. Open MATLAB and navigate to the MoSyn directory.

3. Run the `mosyn` function to start the application:
```matlab
mosyn
```

### Usage

1. Use the GUI to configure your project, including input files, settings, and graph measures.
2. Run the analysis by clicking the corresponding buttons in the GUI.
3. View and export the calculated graph features and measures to output files.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The development team and contributors
- External libraries and resources used in the project
