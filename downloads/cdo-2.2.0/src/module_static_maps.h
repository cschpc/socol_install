#ifndef MODULE_STATIC_MAPS_H
#define MODULE_STATIC_MAPS_H

#include <map>
#include <string>

struct module_t;

/*
 * This file serves no purpose but to make these "hidden" static definitions more visible for
 * future developers.
 */

/**
 * ==========================================================================================
 * Functions Containing static structures for static initaializaton of modules,
 * aliases and operators.  All static maps contained within the following
 * functions are necessary to not encounter problems with the non guaranteed
 * order of initalization of static variables.
 * ==========================================================================================
 */
std::map<std::string, module_t> &get_modules();

std::map<std::string, std::string> &get_aliases();

std::map<std::string, std::string> &get_module_map();
//==========================================================================================
//
#endif
