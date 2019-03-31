//
//  main.h
//  informatiqueGraphique
//
//  Created by amr miguil on 27/03/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//


#include "std_image_write.h"
#include <random>
#include <vector>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <iostream>
#include "math.h"
#include "stb_image.h"
/**
 * Définition d'un moteur de génération de nombre aléatoires mis ici car utilisé par classes.h lors de l'implémentation du comportement des objets de TypeMateriel "Fresnel"
 */
std::vector<std::default_random_engine> listEngines;
static  std::default_random_engine engine;
static  std::uniform_real_distribution<double> distrib(0, 1);

#include "classes.h"
#include "SceneObject.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Geometry.h"
#include "Scene.h"
#include <omp.h>
#include <chrono>
