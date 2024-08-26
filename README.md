# STATCOM in Multi-Area Multi-Machine System to Improve Low-Frequency Control and Voltage Stability

This repository contains the implementation of a Hybrid Artificial Bee Colony and Particle Swarm Optimization (HABC-PSO) approach for Static Synchronous Compensator (STATCOM) control in a multi-area, multi-machine power system. The goal is to improve system damping, address low-frequency and sub-synchronous oscillations, and enhance overall power system stability.

## Project Description

The proposed control strategy utilizes a hybrid optimization technique, combining Artificial Bee Colony (ABC) and Particle Swarm Optimization (PSO), to fine-tune the proportional-integral (PI) controllers of the STATCOM's AC and DC voltage regulators. The optimization is performed using a time-domain objective function to ensure optimal performance under various operating conditions.

### Key Features

- **HABC-PSO Optimization**: A novel hybrid approach integrating ABC and PSO algorithms to optimize the STATCOM control parameters.
- **Improved Damping and Stability**: The optimized STATCOM enhances system damping, reducing low-frequency and sub-synchronous oscillations.
- **Comprehensive Analysis**: Both steady-state and transient analyses are conducted using a Kundur multi-area four-machine model.
- **Performance Evaluation**: The simulation results include assessments of generator rotor speed deviation, angle deviation, and bus-bar voltage deviation. The effectiveness of the STATCOM is evaluated by comparing scenarios with and without the STATCOM.

### Simulation Results

The simulation outcomes demonstrate that the HABC-PSO-based STATCOM effectively mitigates inter-area oscillations, providing a more resilient and efficient solution for maintaining power system stability. The results show reduced settling times for rotor speed deviation (Δω) and angle deviation (Δδ) of generators across different areas, affirming the robustness of the proposed control scheme.

## Tools and Technologies

- MATLAB Simulink
- Power System Toolbox
- Hybrid Optimization Algorithms (HABC-PSO)

## Repository Contents

- **/models**: Contains the Simulink models used for the simulations.
- **/scripts**: MATLAB scripts for running simulations and analyzing results.
- **/results**: Output data and plots from the simulations.
- **README.md**: This file.

## Getting Started

To get started with the project, clone the repository and open the Simulink models in MATLAB. Follow the instructions in the `/scripts` directory to run the simulations and analyze the results.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or collaboration, please contact [Your Name](mailto:your-email@example.com).
