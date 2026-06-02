# Contributors Guide

Thank you for considering contributing to WAVI.jl! If you're interested in contributing to the development of WAVI we want your help no matter how big or small a contribution you make! 

## Table of Contents
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Code Contribution](#code-contribution)
- [Pull Request Process](#pull-request-process)
- [License](#license)

## How to Contribute

### Reporting Bugs

If you encounter a bug, please help us fix it by following these steps:

1. Ensure the bug is not already reported by checking the [issue tracker](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues).
2. If the bug isn't reported, open a new issue. Clearly describe the issue, including steps to reproduce it.

### Suggesting Enhancements

If you have ideas for enhancements, new features, or improvements, we'd love to hear them! Follow these steps:

1. Check the [issue tracker](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues) to see if your suggestion has been discussed.
2. If not, open a new issue, providing a detailed description of your suggestion and the use case it addresses.

### Code Contribution

If you'd like to contribute code to the project:

1. Fork the repository.
2. Clone your fork: `git clone https://github.com/WAVI-ice-sheet-model/WAVI.jl`
3. Create a new branch for your changes: `git checkout -b feature-branch`
4. Make your changes and commit them with a clear message.
5. Push your changes to your fork: `git push origin feature-branch`
6. Open a pull request against the `main` branch of the main repository.


## Local development

To set up WAVI for local development (e.g. to contribute or test changes), follow these steps:

**Clone the repository and create branch as needed (above section)**:

```bash
git clone https://github.com/WAVI-ice-sheet-model/WAVI.jl.git
cd WAVI.jl
```

Lets say this is cloned to: `/git/WAVI.jl`.

### Editable install

If you would like to use your local WAVI code in another project (i.e. outside of the WAVI repo), do the following from that project directory:

```julia
julia> ]

(@v1.12) pkg>activate .

(@v1.12) pkg>develop /git/WAVI.jl
```

This adds WAVI to your current project environment in editable mode, pointing to your local clone. Any changes you make to the code in `/git/WAVI.jl` will take effect immediately — just restart your Julia session or re-include the relevant files to see the updates.

Make sure to update `/git/WAVI.jl` in the above example to match the actual location where you cloned the `WAVI.jl` repo.

## Pull Request Process

Please ensure your pull request follows these guidelines:

1. Adheres to the coding standards.
2. Includes relevant tests for new functionality.
3. Has a clear commit history and messages.
4. References the relevant issue if applicable.

Please don't hesistate to get in touch to discuss, or with any questions!

## License

By contributing to this project, you agree that your contributions will be licensed under the [LICENSE](https://github.com/WAVI-ice-sheet-model/WAVI.jl/blob/main/LICENSE) file of this repository.

