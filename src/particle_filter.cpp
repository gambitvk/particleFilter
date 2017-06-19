/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <utility>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    default_random_engine gen;
    num_particles = 500;
    normal_distribution<double> dist_x(x,std[0]);
    normal_distribution<double> dist_y(y,std[1]);
    normal_distribution<double> dist_psi(theta,std[2]);

    Particle tmp;

    for(int i = 0; i < num_particles ; i++)
    {
        tmp.id  = i;
        tmp.x = dist_x(gen);
        tmp.y = dist_y(gen);
        tmp.theta = dis_psi(gen);
        tmp.weight = 1.0;
        particles.push_back(tmp);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    default_random_engine gen;
    normal_distribution<double> dist_x(0,std[0]);
    normal_distribution<double> dist_y(0,std[1]);
    normal_distribution<double> dist_psi(0,std[2]);

    for (int i = 0; i < num_particles; i++) 
    {
        if (fabs(yaw_rate) < 0.001) 
        {  
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        } 
        else 
        {
            particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
            particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
            particles[i].theta += yaw_rate * delta_t;
        }

        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_psi(gen);
    }
}
    
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    //
}

//the fact that I can only submit particle_filter.cpp means I can not change any of the existing structure nor I could add helper functions as part of the class.. so I can only add a non class function here and use it in my class fucntion above. 
//
//PLEASE do not complain about my c++ style here as I had no choice. 
//

void transform( const LandmarkObs & obs, const Particle & par, double &t_x, double &t_y)
{
    t_x = cos(par.theta) * obs.x - sin(par.theta) * obs.y + par.x;
    t_y = sin(par.theta) * obs.x + cos(par.theta) * obs.y + par.y;
}

void updateOneParticle(std::vector<LandmarkObs> &predicted, std::vector<LandmarkObs> & observations
                   ,Particle &par, double &std_x, double &std_y)
    std::pair<int,double> min{-1,0.0};
    std::pair<double,double> xy{0.0,0.0};
    std::pair<double,double> p_xy{0.0,0.0};
    double dis = 0.0;
    par.weight = 1.0;
    for(const auto &obs : observations)
    {
        for(const auto &pred : predicted)
        {
            transform(obs,par,xy.first,xy.second);
            dis = dist(xy.first,xy.second,pred.x,pred.y);
            if (min.first == -1 || dis < min.second) 
            {
                min.first = pred.id;
                min.second = dis;
                p_xy.first = pred.x;
                p_xy.second = pred.y;
            }
        }

        //update weight for particle here
        par.weight *=  exp( -( (pow(p_xy.first-xy.first,2)/(2*pow(std_x, 2))) 
                              +(pow(p_xy.second-xy.second,2)/(2*pow(std_y, 2)))
                             )
                           ) / (2*M_PI*std_x*std_y) ;

        min.first = -1;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks)
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    weights.clear();
    for(auto &par : particles)
    {
        updateOneParticle(map_landmarks.landmark_list,observations,par
                          ,std_landmark[0],std_landmark[1] ); 
        weights.push_back(par.weight);
    }
    
}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::discrete_distribution<> dist_w(weights.begin(), weights.end());
  default_random_engine gen;

 for (int p = 0; p<num_particles; p++) 
 {
  dist_w(gen);
  particles[densityVector(gen)].weights = -1.0;
 }

 auto notPicked = [](auto &a) {return a.weights != -1.0 ;};
 std::remove_if( particles.begin(),particles.end(),notPicked);
 num_particles = particles.size();
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

