/****************************************************************************
 * Copyright ©2017 Brian Curless.  All rights reserved.  Permission is hereby
 * granted to students registered for University of Washington CSE 457 or CSE
 * 557 for use solely during Autumn Quarter 2017 for purposes of the course.
 * No other use, copying, distribution, or modification is permitted without
 * prior written consent. Copyrights for third-party components of this work
 * must be honored.  Instructors interested in reusing these course materials
 * should contact the author.
 ****************************************************************************/
#include "particlesystem.h"
#include "spherecollider.h"
#include "planecollider.h"
#include "cylindercollider.h"
#include <scene/sceneobject.h>

REGISTER_COMPONENT(ParticleSystem, ParticleSystem)

ParticleSystem::ParticleSystem() :
    ParticleGeometry({"Sphere"}, 0),
    ParticleMaterial(AssetType::Material),
    Mass(0.1f, 0.0f, 10.0f, 0.1f),
    Period(0.5f, 0.0f, 1.0f, 0.01f),
    InitialVelocity(glm::vec3(5.0f, 5.0f, 0.0f)),
    ConstantF(glm::vec3(0.0f, -9.8f, 0.0f)),
    DragF(0.0f, 0.0f, 10.0f, 0.01f),
    constant_force_(ConstantF.Get()),
    // REQUIREMENT:
    // init drag force with DragF -- refer to how we deal with constant_force_
    // remember  (f = -k_d * v), where DragF represents k_d
    drag_force_(DragF.Get()),
    num_particles_(0),
    particle_index_(0),
    simulating_(false)
{
    AddProperty("Geometry", &ParticleGeometry);
    AddProperty("Material", &ParticleMaterial);
    AddProperty("Mass", &Mass);
    AddProperty("Period (s)", &Period);
    AddProperty("Initial Velocity", &InitialVelocity);
    AddProperty("Constant Force", &ConstantF);
    AddProperty("Drag Coefficient", &DragF);

    ParticleGeometry.ValueSet.Connect(this, &ParticleSystem::OnGeometrySet);

    forces_.push_back(std::shared_ptr<Force>(&constant_force_));

    // REQUIREMENT: 
    //    add viscous drag force into forces_ array if your drag force also inherit from class Force
    //    If not, you could use your own way to prepare your drag force
    forces_.push_back(std::shared_ptr<Force>(&drag_force_));
}

void ParticleSystem::UpdateModelMatrix(glm::mat4 model_matrix) {
   model_matrix_ = model_matrix;
}

void ParticleSystem::EmitParticles() {
    if (!simulating_) return;

    // REQUIREMENT:
    // Create some particles!
    //    - We have designed a class Particle for you
    //    - We've provided some UI controls for you
    //          -- Mass.Get() defines particle mass, and
    //          -- InitialVelocity.Get() defines particle init velocity in local object space
    //    - Notice particles should be created in world space. (use model_matrix_ to transform particles from local object space to world space)
    // Store particles in the member variable particles_
    // For performance reasons, limit the amount of particles that exist at the same time
    // to some finite amount (MAX_PARTICLES). Either delete or recycle old particles as needed.
    glm::vec3 init_vel_world = (glm::mat3(model_matrix_)*InitialVelocity.Get());
    glm::vec3 init_pos_world = (model_matrix_*glm::vec4(0.0f, 0.0f, 0.0f, 1.0f)).xyz;
    glm::vec3 init_rot_world = (model_matrix_*glm::vec4(0.0f, 0.0f, 0.0f, 1.0f)).xyz;
    Particle* p = new Particle(Mass.Get(), init_pos_world, init_vel_world, init_rot_world);
    particles_.push_back(std::unique_ptr<Particle>(p));
    if (particles_.size() > MAX_PARTICLES) {
        particles_.erase(particles_.begin());
    }
    // Reset the time
    time_to_emit_ = Period.Get();
}

std::vector<Particle*> ParticleSystem::GetParticles() {
    // Return a vector of particles (used by renderer to draw them)
    std::vector<Particle*> particles;
    for (auto& particle : particles_) particles.push_back(particle.get());
    return particles;
}

void ParticleSystem::StartSimulation() {
    simulating_ = true;
    constant_force_.SetForce(ConstantF.Get());
    // REQUIREMENT:
    // Set your added drag force as DragF.Get() -- Refer to what we did on constact_force_
    drag_force_.SetForce(DragF.Get());
    ResetSimulation();
}

void ParticleSystem::UpdateSimulation(float delta_t, const std::vector<std::pair<SceneObject*, glm::mat4>>& colliders) {
    if (!simulating_) return;

    // Emit Particles
    time_to_emit_ -= delta_t;
    if (time_to_emit_ <= 0.0) EmitParticles();

   // REQUIREMENT: 
   // For each particle ...
   //      Calculate forces
   //      Solve the system of forces using Euler's method
   //      Update the particle
   //      Check for and handle collisions
    for (auto& p : particles_) {
        glm::vec3 force_accum = {0, 0, 0};
        for (auto& f : forces_) {
            force_accum += f->GetForce(*p);
        }
        p->Position += delta_t * p->Velocity;
        p->Velocity += delta_t * force_accum/p->Mass;
    }

   // Collision code might look something like this:
   for (auto& kv : colliders) {
       SceneObject* collider_object = kv.first;
       glm::mat4 collider_model_matrix = kv.second;
       glm::mat4 inv_collider_model_matrix = inverse(collider_model_matrix);

       static const double EPSILON = 0.1;
       float particle_radius = 0.5f;

       for (auto& p : particles_) {
           glm::vec3 p_pos = (inv_collider_model_matrix*glm::vec4(p->Position.x, p->Position.y, p->Position.z, 1.0f)).xyz;
           glm::vec3 p_vel = glm::mat3(inv_collider_model_matrix)*p->Velocity;
           // When checking collisions, remember to bring particles from world space to collider local object space
           // The trasformation matrix can be derived by taking invese of collider_model_matrix
           if (SphereCollider* sphere_collider = collider_object->GetComponent<SphereCollider>()) {
               // Check for Sphere Collision
               glm::vec3 N = p_pos/length(p_pos);
               if (length(p_pos) < sphere_collider->Radius.Get() + particle_radius + EPSILON && dot(p_vel, N) < 0) {
                   glm::vec3 vel_norm = dot(N, p_vel)*N;
                   glm::vec3 vel_tan = p_vel - vel_norm;
                   p_vel = vel_tan - glm::vec3(sphere_collider->Restitution.Get()*vel_norm.x, sphere_collider->Restitution.Get()*vel_norm.y, sphere_collider->Restitution.Get()*vel_norm.z);
                   p->Velocity = glm::mat3(collider_model_matrix)*p_vel;
               }
           } else if (PlaneCollider* plane_collider = collider_object->GetComponent<PlaneCollider>()) {
               // Check for Plane Collision
               glm::vec3 N = {0, 0, 1.0f};
               if (dot(p_pos, N) <= particle_radius + EPSILON && dot(p_pos, N) > particle_radius - EPSILON && dot(p_vel, N) < 0) {
                   if (p_pos.x > -plane_collider->Width.Get()/2.0f && p_pos.x < plane_collider->Width.Get()/2.0f &&
                       p_pos.y > -plane_collider->Height.Get()/2.0f && p_pos.y < plane_collider->Height.Get()/2.0f) {
                       glm::vec3 vel_norm = dot(N, p_vel)*N;
                       glm::vec3 vel_tan = p_vel - vel_norm;
                       p_vel = vel_tan - glm::vec3(plane_collider->Restitution.Get()*vel_norm.x, plane_collider->Restitution.Get()*vel_norm.y, plane_collider->Restitution.Get()*vel_norm.z);
                       p->Velocity = glm::mat3(collider_model_matrix)*p_vel;
                   }
               }
           }
           // When updating particle velocity, remember it's in the worls space.
       }
   }
}

void ParticleSystem::StopSimulation() {
    simulating_ = false;
}

void ParticleSystem::ResetSimulation() {
    // Clear all particles
    particles_.clear();
    time_to_emit_ = Period.Get();
}

bool ParticleSystem::IsSimulating() {
    return simulating_;
}


void ParticleSystem::OnGeometrySet(int c) {
    GeomChanged.Emit(ParticleGeometry.GetChoices()[c]);
}
