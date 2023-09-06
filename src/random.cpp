
#include "header.hpp"
#include "model.hpp"

using namespace std;

double Model::random_real(double min, double max)
{
  return uniform_real_distribution<>(min, max)(gen);
}

double Model::random_normal(double sigma)
{
  return normal_distribution<>(0., sigma)(gen);
}

unsigned Model::random_geometric(double p)
{
  return geometric_distribution<>(p)(gen);
}

unsigned Model::random_unsigned()
{
  return gen();
}

void Model::InitializeRandomNumbers()
{
  if(not set_seed)
  {
    // 'Truly random' device to generate seed
    std::random_device rd;
    seed = rd();
  }

  gen.seed(seed);
}

unsigned Model::random_poisson(double lambda)
{
  return poisson_distribution<>(lambda)(gen);
}

int Model::random_exponential(double lambda)
{
  return exponential_distribution<>(lambda)(gen);
}
