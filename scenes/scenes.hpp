// Scenes shared by the tests, the benchmark and the demo.
#pragma once

#include <string>
#include <vector>

#include "rbe/world.hpp"

namespace rbe::scenes {

BodyId addGround(World& w);
// Tower of n unit boxes resting on the ground. Returns the bodies bottom to top.
std::vector<BodyId> stack(World& w, int n);
// Pyramid of unit boxes, `rows` boxes wide at the base.
std::vector<BodyId> pyramid(World& w, int rows);
// n unit boxes dropped into a walled pit as a loose grid so they pile up.
std::vector<BodyId> pile(World& w, int n, unsigned seed = 1);
// One box with three different moments of inertia, gravity off, spinning fast about its middle axis
// (the axis that flips over: tennis racket / Dzhanibekov effect).
BodyId spinningBox(World& w, float spinRadPerSec = 20.0f);
std::vector<BodyId> dominoes(World& w, int n);
std::vector<BodyId> sphereRain(World& w, int n, unsigned seed = 1);

// Builds a scene by name ("stack", "pyramid", "pile", "spin", "dominoes", "spheres"). Returns false if unknown.
bool build(World& w, const std::string& name, int count);

}  // namespace rbe::scenes
