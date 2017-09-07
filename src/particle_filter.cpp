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
#include <assert.h>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;

	num_particles = 50;

	//Based on x, y and theta generate initial particles with standard deviation std[0], std[1] and std[2]
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	for (int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = x + dist_x(gen);
		particle.y = y + dist_y(gen);
		particle.theta = theta + dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	double x, xf = 0.0;
	double y, yf = 0.0;
	double theta, thetaf = 0.0;

	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	default_random_engine gen;

	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);


	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles.at(i);

		x = particle.x;
		y = particle.y;
		theta = particle.theta;

		if (fabs(yaw_rate) > 0.00001)
		{
			xf = x + (velocity/yaw_rate)*( sin(theta + yaw_rate * delta_t) - sin(theta));
			yf = y + (velocity/yaw_rate)* (cos(theta) - cos(theta + yaw_rate * delta_t));
			thetaf = theta + yaw_rate * delta_t;
		}
		else
		{
			xf = x + velocity * delta_t * cos(theta);
			yf = y + velocity * delta_t * sin(theta);
		}

		//Add noise
		particle.x = xf + dist_x(gen);
		particle.y = yf + dist_y(gen);

		particle.theta = thetaf + dist_theta(gen);
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < observations.size(); i++)
	{
		double minDist = numeric_limits<double>::max();

		for (int j= 0; j < predicted.size(); j++)
		{
			double distance = dist(	observations[i].x, observations[i].y,
									predicted[j].x, predicted[j].y);

			if (distance < minDist)
			{
				observations[i].id = predicted[j].id;
				minDist = distance;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles.at(i);

		vector<LandmarkObs> validLandmarks;

		//Iterate through all the landmarks
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++)
		{
			const Map::single_landmark_s& landMark = map_landmarks.landmark_list[k];

			//if ((fabs(particle.x - landMark.x_f) <= sensor_range) &&
			//		(fabs(particle.y - landMark.y_f) <= sensor_range))

			if (dist(particle.x, particle.y, landMark.x_f, landMark.y_f) <= sensor_range)
			{

				validLandmarks.push_back(LandmarkObs{landMark.id_i, landMark.x_f, landMark.y_f});
			}
		}


		vector<LandmarkObs> transObs;

		for (int k = 0; k < observations.size(); k++)
		{
		      double tx = cos(particle.theta)*observations[k].x - sin(particle.theta)*observations[k].y + particle.x;
		      double ty = sin(particle.theta)*observations[k].x + cos(particle.theta)*observations[k].y + particle.y;
		      transObs.push_back(LandmarkObs{observations[k].id, tx, ty});
		}

		dataAssociation(validLandmarks, transObs);

		particle.weight = 1.0;

		for (int k = 0; k < transObs.size(); k++)
		{
			double x = transObs[k].x;
			double y = transObs[k].y;
			double ux = 0.0;
			double uy = 0.0;
			double stdx = std_landmark[0];
			double stdy = std_landmark[1];

			//Get associate landMark
			for (int j = 0; j < validLandmarks.size();j++)
			{
				if (transObs[k].id == validLandmarks[j].id)
				{
					ux = validLandmarks[j].x;
					uy = validLandmarks[j].y;
					break;
				}
			}

			double obsW = ( 1/(2* M_PI * stdx * stdy)) * exp( -( pow(x-ux,2)/(2*pow(stdx, 2)) + (pow(y-uy,2)/(2*pow(stdy, 2)))));

			particle.weight *= obsW;
		}

		weights.push_back(particle.weight);

		cout << "Particle Id: "<< particle.id <<" Weight: " << particle.weight << endl;

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


    std::random_device rd;
    std::mt19937 gen(rd());

    std::discrete_distribution<> d(weights.begin(), weights.end());

    std::vector<Particle> resampled_particle_list;


    for (int i = 0; i < particles.size(); i++)
    {

    	int w = d(gen);

    	cout << "w: " << w;

    	Particle part(particles[w]);

    	Particle& p = particles.at(w);

    	cout << " part id: " << part.id << " p id:"<< p.id << endl;

    	resampled_particle_list.push_back(part);
    }

    particles.clear();

    particles = resampled_particle_list;

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


void ParticleFilter::getMapCoordinates(const Particle& particle, const double xc, const double yc, double& xm, double& ym)
{
	xm = particle.x + (cos (particle.theta)* xc) - (sin (particle.theta)* yc);
	ym = particle.y + (sin (particle.theta)* xc) + (cos (particle.theta)* yc);
}

LandmarkObs ParticleFilter::getNearestObservationForLandmark(const Particle& particle, const Map::single_landmark_s& landMark, const std::vector<LandmarkObs>& landmarkObs)
{
	double 		minDistance = 0.0;
	LandmarkObs selectedObs = landmarkObs[0];


	assert (landmarkObs.size() != 0);

	for (int j = 0; j < landmarkObs.size(); j++)
	{
		//Get map coordinates for this observation for this particle
		double xm = 0.0;
		double ym = 0.0;
		double distance = 0.0;

		getMapCoordinates(particle, landmarkObs[j].x, landmarkObs[j].y, xm, ym);

		distance = dist(xm, ym, landMark.x_f, landMark.y_f);

		//cout << "minDistance: " << minDistance << " distance: " << distance <<endl;

		if (j == 0)
		{
			minDistance = distance;
			selectedObs.id = landmarkObs[j].id;
			selectedObs.x = xm;
			selectedObs.y = ym;
			continue;
		}

		if (minDistance > distance)
		{
			minDistance = distance;
			selectedObs.id = landmarkObs[j].id;
			selectedObs.x = xm;
			selectedObs.y = ym;
		}
	}

	return selectedObs;
}


const Map::single_landmark_s ParticleFilter::get_landmark_from_id(int id, const Map& mapLandmarks)
{
	Map::single_landmark_s landMark;

	for (int i = 0; i < mapLandmarks.landmark_list.size(); i++)
	{
		if (mapLandmarks.landmark_list[i].id_i == id)
		{
			Map::single_landmark_s lMark = {mapLandmarks.landmark_list[i].id_i,
											mapLandmarks.landmark_list[i].x_f,
											mapLandmarks.landmark_list[i].y_f};
			return lMark;
		}
	}

	assert (false);

	return landMark;
}


Particle ParticleFilter::getParticle(int particleId)
{

	//ToDo function should be modified to return error if particleId is not found
	Particle particle;
	bool bFound =  false;

	for (int i = 0; i < particles.size(); i++)
	{
		if (particles[i].id == particleId)
		{
			return Particle(particles[i]);
		}
	}

	assert (bFound);

	return particle;
}
