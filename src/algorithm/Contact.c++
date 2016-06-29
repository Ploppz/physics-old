/* src */
#include "Contact.h"
#include <debug.h>
#include <geometry/Polygon.h>
#include <geometry/geometry.h>

DepthContact simple_contact_conversion(TimeContact& in_contact)
{
    DepthContact out;

    // Identify edge & vertex of the contact
    EdgePoint& edge_point = in_contact.ref_point;
    EdgePoint& vertex_point = in_contact.subj_point;
    if (vertex_point.alpha != 0)
        std::swap(edge_point, vertex_point);
    //
    Polygon::Edge edge(edge_point.index, edge_point.parent);
    out.depth = distance_line(vertex_point.point_t(), edge.start_tr(), edge.end_tr());
    out.normal = in_contact.normal;
    out.subj_point = vertex_point;
    out.ref_point = edge_point;
    return out;
}
void DepthContact::ensure_normal_toward_subject()
{
    DebugBegin();
    float d = glm::dot(subj_point.point_t() - ref_point.point_t(), normal);
    // dout << "dot = " << d / (glm::length(subj_point.point_t() - ref_point.point_t()) * glm::length(normal)) << newl;
    if (d < 0) {
        normal = - normal;
        // dout << "Flipping normal" << newl;
    }
}
