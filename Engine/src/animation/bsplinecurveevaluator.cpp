/****************************************************************************
 * Copyright ©2017 Brian Curless.  All rights reserved.  Permission is hereby
 * granted to students registered for University of Washington CSE 457 or CSE
 * 557 for use solely during Autumn Quarter 2017 for purposes of the course.
 * No other use, copying, distribution, or modification is permitted without
 * prior written consent. Copyrights for third-party components of this work
 * must be honored.  Instructors interested in reusing these course materials
 * should contact the author.
 ****************************************************************************/
#include "bsplinecurveevaluator.h"

std::vector<glm::vec2> BSplineCurveEvaluator::EvaluateCurve(const std::vector<glm::vec2> &ctrl_pts, int density) const {
    std::vector<glm::vec2> evaluated_pts;

    // REQUIREMENT:
    // Currently this function returns points for a Linear Evaluation.
    // Replace this code with code that returns evaluated points for a B-Spline
    // curve. Be sure to respect the extend_x_ and wrap_ flags.

    if (density == 0) density = 100;
    if (ctrl_pts.size() < 3) {
        for (size_t i = 0; i < ctrl_pts.size()-1; i++) {
            for (int j = 0; j < density; j++) {
                float t = j/(float) density;
                glm::vec2 p = t*ctrl_pts[i+1] + (1-t)*ctrl_pts[i];
                evaluated_pts.push_back(p);
            }
        }
        evaluated_pts.push_back(ctrl_pts.back());
    } else {
        for (int j = 0; j < density; j++) {
            float t = j/(float) density;
            glm::vec2 p = 1.0f/6.0f*((-t*t*t + 3*t*t - 3*t + 1)*ctrl_pts[0] + (3*t*t*t - 6*t*t + 4)*ctrl_pts[0] + (-3*t*t*t + 3*t*t + 3*t + 1)*ctrl_pts[0] + (t*t*t)*ctrl_pts[1]);
            evaluated_pts.push_back(p);
        }
        for (int j = 0; j < density; j++) {
            float t = j/(float) density;
            glm::vec2 p = 1.0f/6.0f*((-t*t*t + 3*t*t - 3*t + 1)*ctrl_pts[0] + (3*t*t*t - 6*t*t + 4)*ctrl_pts[0] + (-3*t*t*t + 3*t*t + 3*t + 1)*ctrl_pts[1] + (t*t*t)*ctrl_pts[2]);
            evaluated_pts.push_back(p);
        }
        for (size_t i = 0; i < ctrl_pts.size()-3; i++) {
            for (int j = 0; j < density; j++) {
                float t = j/(float) density;
                glm::vec2 p = 1.0f/6.0f*((-t*t*t + 3*t*t - 3*t + 1)*ctrl_pts[i] + (3*t*t*t - 6*t*t + 4)*ctrl_pts[i+1] + (-3*t*t*t + 3*t*t + 3*t + 1)*ctrl_pts[i+2] + (t*t*t)*ctrl_pts[i+3]);
                evaluated_pts.push_back(p);
            }
        }
        for (int j = 0; j < density; j++) {
            float t = j/(float) density;
            glm::vec2 p = 1.0f/6.0f*((-t*t*t + 3*t*t - 3*t + 1)*ctrl_pts[ctrl_pts.size()-3] + (3*t*t*t - 6*t*t + 4)*ctrl_pts[ctrl_pts.size()-2] + (-3*t*t*t + 3*t*t + 3*t + 1)*ctrl_pts[ctrl_pts.size()-1] + (t*t*t)*ctrl_pts[ctrl_pts.size()-1]);
            evaluated_pts.push_back(p);
        }
        for (int j = 0; j < density; j++) {
            float t = j/(float) density;
            glm::vec2 p = 1.0f/6.0f*((-t*t*t + 3*t*t - 3*t + 1)*ctrl_pts[ctrl_pts.size()-2] + (3*t*t*t - 6*t*t + 4)*ctrl_pts[ctrl_pts.size()-1] + (-3*t*t*t + 3*t*t + 3*t + 1)*ctrl_pts[ctrl_pts.size()-1] + (t*t*t)*ctrl_pts[ctrl_pts.size()-1]);
            evaluated_pts.push_back(p);
        }
    }
    if (extend_x_) ExtendX(evaluated_pts, ctrl_pts);
    return evaluated_pts;
}
