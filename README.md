# Towards Joint Neutrino–Seismic Tomography

## Introduction

The Earth's deep interior still holds many mysteries. Seismic imaging has revealed its structure in increasing detail, but seismic velocities alone cannot tell us whether a given anomaly is thermal, chemical, or both. To help resolve this ambiguity, we explore **atmospheric neutrinos** as a complementary **probe of the deep Earth**. 

**Neutrinos** are nearly massless particles that interact extremely weakly with matter. They come in three *flavours* (electron, muon and tau), each with a corresponding antiparticle (antineutrino). Atmospheric neutrinos are produced continuously by cosmic-ray interactions in the upper atmosphere and cross the entire planet almost unimpeded. On the way, they undergo **flavour oscillations**: a neutrino created with one flavour has a probability of being detected as another, and this probability depends on the electron density of the matter along its path.

Because **electron density** depends on both mass density and chemical composition, neutrino observations offer a sensitivity to the deep Earth that is fundamentally different from seismic waves. 

Our paper, [Sensitivity of neutrino oscillations to the Earth’s interior properties](https://www.sciencedirect.com/science/article/pii/S0012821X26005340) (Isabel A. Goos, Nobuaki Fuji, Stéphanie Durand, Yael Deniz Véronique Van Elewyck, João A. B. Coelho, Eric Mittelstaedt), establishes the theoretical framework needed to include neutrino oscillation measurements in a future joint inversion with seismic data. 

## Repository Overview
 
This repository gathers subsequent code and studies I developed for this purpose. Each folder has its own `README.md` with more detail.
 
| Folder | Contents |
|---|---|
| [`NeutrinoInversion/`](NeutrinoInversion/) | Inversion approaches using neutrino data alone, with successively more advanced models, and tests on synthetic neutrino data. |
| [`ForwardModel/`](ForwardModel/) | Core components of the neutrino forward model. They are coded and tested individually here before being integrated into the architecture of the [Flexopt](https://github.com/seismobassoon/flexOPT.jl) framework. |
| [`BinsStudy/`](BinsStudy/) | Study of the binning scheme for neutrino data. |

## Collaboration
 
This is an interdisciplinary effort between:
 
- **Institut de Physique du Globe de Paris (IPGP):** Nobuaki Fuji and myself
- **AstroParticule et Cosmologie (APC), Paris:** Véronique Van Elewyck and João A. B. Coelho
- **Laboratoire de Géologie de Lyon (LGL-TPE):** Stéphanie Durand
- **University of Idaho:** Yael Deniz and Eric Mittelstaedt
