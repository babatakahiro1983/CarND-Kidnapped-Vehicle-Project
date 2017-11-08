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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	/**
	* init Initializes particle filter by initializing particles to Gaussian
	*   distribution around first position and all the weights to 1.
	* @param x Initial x position [m] (simulated estimate from GPS)
	* @param y Initial y position [m]
	* @param theta Initial orientation [rad]
	* @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	*   standard deviation of yaw [rad]]
	*/

	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	
	// Set standard deviations for x, y, and theta
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	// Creates a normal (Gaussian) distribution for x, y and theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);


	for (int loop_1 = 0; loop_1 < num_particles; ++loop_1) {

		Particle p_tmp;
		p_tmp.id = loop_1;
		p_tmp.x = dist_x(gen);
	        p_tmp.y = dist_y(gen);
		p_tmp.theta = dist_theta(gen);
		p_tmp.weight = 1.0;
                particles.push_back(p_tmp);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	/**
	* prediction Predicts the state for the next time step
	*   using the process model.
	* @param delta_t Time between time step t and t+1 in measurements [s]
	* @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	*   standard deviation of yaw [rad]]
	* @param velocity Velocity of car from t to t+1 [m/s]
	* @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	*/

	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];

	for (int loop_1= 0; loop_1 < num_particles; ++loop_1) {

		// Creates a normal (Gaussian) distribution for 
		normal_distribution<double> dist_x(particles[loop_1].x, std_x);
		normal_distribution<double> dist_y(particles[loop_1].y, std_y);
		normal_distribution<double> dist_theta(particles[loop_1].theta, std_theta);

		particles[loop_1].x = dist_x(gen) + velocity / yaw_rate * (sin(dist_theta(gen) + yaw_rate * delta_t) - sin(dist_theta(gen)));
		particles[loop_1].y = dist_y(gen) + velocity / yaw_rate * (cos(dist_theta(gen)) - cos(dist_theta(gen) + yaw_rate * delta_t));
		particles[loop_1].theta = dist_theta(gen) + yaw_rate * delta_t;
	}
}

void ParticleFilter::dataAssociation(const std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	/**
	* dataAssociation Finds which observations correspond to which landmarks (likely by using
	*   a nearest-neighbors data association).
	* @param predicted Vector of predicted landmark observations
	* @param observations Vector of landmark observations
	*/

	double dist_tmp, dist_min;
	dist_min = 100;

	// 観測点ごとにランドマークの最近点を探索する。
	for (int loop_3 = 0; loop_3 < observations.size(); ++loop_3) {
		for (int loop_4 = 0; loop_4 < predicted.size(); ++loop_4) {
			
			dist_tmp = dist(predicted[loop_4].x, predicted[loop_4].y, observations[loop_3].x, observations[loop_3].y);

			if (dist_tmp < dist_min) {
		
				observations[loop_3].id = loop_4; 
				dist_min = dist_tmp;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	double sig_x, sig_y, gauss_norm, exponent, mu_x, mu_y, weight;
	std::vector<LandmarkObs> observations_map, landmarks;

	
	for (int loop_1 = 0; loop_1 < num_particles; ++loop_1) {

		particles[loop_1].weight = 1.0;		

		// (1)Transform
		for (int loop_2 = 0; loop_2 < observations.size(); ++loop_2) {

			LandmarkObs obs_tmp;

			obs_tmp.id = loop_2;
			
			// transform to map x coordinate
			obs_tmp.x = particles[loop_1].x + (cos(particles[loop_1].theta) * observations[loop_2].x) - (sin(particles[loop_1].theta) * observations[loop_2].y);

			//transform to map y coordinate
			obs_tmp.y = particles[loop_1].y + (sin(particles[loop_1].theta) * observations[loop_2].x) + (cos(particles[loop_1].theta) * observations[loop_2].y);

			observations_map.push_back(obs_tmp);

		}
		cout << "landmarks " << observations_map.size() << endl;

		// (2)Associate
		for (int loop_4 = 0; loop_4 < map_landmarks.landmark_list.size(); ++loop_4) {

			LandmarkObs landmarks_tmp;

			landmarks_tmp.id = map_landmarks.landmark_list[loop_4].id_i;
			landmarks_tmp.x = map_landmarks.landmark_list[loop_4].x_f;
			landmarks_tmp.y = map_landmarks.landmark_list[loop_4].y_f;

			landmarks.push_back(landmarks_tmp);
			
		}
		cout << "landmarks " << landmarks.size() << endl;

		dataAssociation(landmarks, observations_map);


		
		
		// (3)updateWeight
		for (int loop_2 = 0; loop_2 < observations.size(); ++loop_2) {
			
			sig_x = std_landmark[0];
			sig_y = std_landmark[1];
			mu_x = map_landmarks.landmark_list[observations[loop_2].id].x_f; // assosiateされたランドマークの座標
			mu_y = map_landmarks.landmark_list[observations[loop_2].id].y_f; // assosiateされたランドマークの座標

			// calculate normalization term
			gauss_norm = (1 / (2 * M_PI * sig_x * sig_y));

			// calculate exponent
			exponent = std::pow((observations_map[loop_2].x - mu_x), 2.0) / std::pow(2 * sig_x, 2.0) 
				     + std::pow((observations_map[loop_2].y - mu_y), 2.0) / std::pow(2 * sig_y, 2.0);

			// calculate weight using normalization terms and exponent
			weight = gauss_norm * exp(-exponent);

			// multiply all the calculated measurement probabilities
			particles[loop_1].weight = particles[loop_1].weight * weight;
		}

		cout << "weight " << particles[loop_1].weight << endl;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rnd;
	std::vector<Particle> particles_temp;
	std::discrete_distribution<int> dist(0,num_particles - 1);
	
	int index = dist(rnd);
	double beta = 0.0;
	double weight_max = 0.0;

	for (int loop_1 = 0; loop_1 < num_particles; ++loop_1) {
		if (weight_max < particles[loop_1].weight) {
			weight_max = particles[loop_1].weight;
		}
	}


	for (int loop_1 = 0; loop_1 < num_particles; ++loop_1) {
			
		Particle p_temp;
		
		beta += rnd() * 2.0 * weight_max;
		while (beta > particles[index].weight) {
			beta -= particles[index].weight;
			index = (index + 1) % num_particles;
		}

		p_temp = particles[index];
		particles_temp.push_back(p_temp);

	}
	
	for (int loop_1 = 0; loop_1 < num_particles; ++loop_1) {

		particles[loop_1] = particles_temp[loop_1];
	}

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
