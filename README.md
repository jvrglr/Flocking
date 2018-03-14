# Flocking
The concept of *emergence* is exemplified through the well-known movement of a flock of birds.
From very simple rules at the individual (bird) level a non trivial movement at the collective (flock) level will give rise.
# Codes
* **Vicsek_model.f**: FORTRAN 77 code, implementation of 2D Vicsek model by using link-cell + ghost cells method. To generate video with evolution of the flock:
 1. Compile and run Vicsek_model.f
 1. run Execute.
 
The output is a .mp4 video with the evolution of the flock starting from a homogeneous initial condition (positions and velocities setted up at random). All frames of the movie will also be saved as .png files, to avoid this, open the Execute file with a text editor and follow the instructions. CAUTION! THE EXECUTE SCRIPT WILL DELETE ALL .DAT FILES!
* **Boid.py**: Python 3 code, implementation of a 2D naive boids model. The output is a .mp4 video with the evolution of the flock starting from a homogeneous initial condition (positions and velocities setted up at random).
## URL addresses
* [Movie0](www.youtube.com/watch?v=xxw3zglK7Os&t=0s&index=16&list=PL7A9POR1j9Mox3jOb0YYVWo4Y_6tCVkeo) www.youtube.com/watch?v=xxw3zglK7Os&t=0s&index=16&list=PL7A9POR1j9Mox3jOb0YYVWo4Y_6tCVkeo - Dr. Andrea Cavagna explains basis of the theory plus important empirical results.
* [Movie1](https://www.youtube.com/watch?v=PTo5Akpjpww&feature=youtu.be) https://www.youtube.com/watch?v=PTo5Akpjpww&feature=youtu.be - Simulated flock of birds by using Boids model.
* Some examples of how a bad parametric choice cand lead to unnatural movements
  * [Movie2](https://www.youtube.com/watch?v=G82fXgspkhs&feature=youtu.be) https://www.youtube.com/watch?v=G82fXgspkhs&feature=youtu.be 
  * [Movie3](https://www.youtube.com/watch?v=UfYhpFwwyQI) https://www.youtube.com/watch?v=UfYhpFwwyQI
  * [Movie4](https://www.youtube.com/watch?v=hFrUGDazobQ&feature=youtu.be) https://www.youtube.com/watch?v=hFrUGDazobQ&feature=youtu.be
  * [Movie5](https://www.youtube.com/watch?v=HG7Wh7oYMis&feature=youtu.be) https://www.youtube.com/watch?v=HG7Wh7oYMis&feature=youtu.be
* [Movie6](https://www.youtube.com/watch?v=KwlSgCWQt4k&feature=youtu.be) https://www.youtube.com/watch?v=KwlSgCWQt4k&feature=youtu.be - Simulated flock of birds by using Vicsek model.
* [Movie7](https://www.youtube.com/watch?v=A8Xt3DlrwhY&feature=youtu.be) https://www.youtube.com/watch?v=A8Xt3DlrwhY&feature=youtu.be - Bats simulated by Boids model. Scene from 1992 Batman Returns.
