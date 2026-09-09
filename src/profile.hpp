// Tracy zones, compiled out unless RBE_TRACY is defined (the demo turns it on).
#pragma once

#ifdef RBE_TRACY
#include <tracy/Tracy.hpp>
#define RBE_ZONE(name) ZoneScopedN(name)
#else
#define RBE_ZONE(name) ((void)0)
#endif
