# ASSA-AVRORA Sieve
A memory-efficient, hybrid CPU-GPU sieve for local prime enumeration in quadratic intervals.

## License
[![License: CC BY-NC 4.0](https://licensebuttons.net/l/by-nc/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc/4.0/)
**Non-commercial use only.**

## Build
```bash
mkdir build && cd build
cmake ..
make
./assa_aurora 67108864
```

## Citation
If you use this code in your research, please cite:
```bibtex
@misc{assa_aurora,
  title={The ASSA-AVRORA Sieve: A Memory-Efficient, Hybrid CPU-GPU Algorithm for Local Prime Enumeration in Quadratic Intervals},
  author={Siarhei Tabalevich, Siarhei Aleksandrov},
  year={2025},
  howpublished={\url{https://github.com/ASSA-NI-ATOM/assa-aurora}},
  note={arXiv preprint}
}
```

## License Details
This work is licensed under **CC BY-NC 4.0** (Creative Commons Attribution-NonCommercial 4.0 International).  
**Commercial use is explicitly prohibited.**  
**Non-commercial research, education, and personal use are freely permitted.**  
All algorithmic components are explicitly placed into the public domain to prevent downstream patenting.  
All commercial rights are reserved by the authors.
