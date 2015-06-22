#include "Scene_polyhedron_selection_item.h"
#include "Scene_polyhedron_selection_item.moc"
struct light_info
{
    //position
    GLfloat position[4];

    //ambient
    GLfloat ambient[4];

    //diffuse
    GLfloat diffuse[4];

    //specular
    GLfloat specular[4];
    GLfloat spec_power;

};

void Scene_polyhedron_selection_item::initialize_buffers(Viewer_interface *viewer)const
{
    //vao containing the data for the unselected facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_facets.data(), positions_facets.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();



        buffers[1].bind();
        buffers[1].allocate(normals.data(), normals.size()*sizeof(float));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[1].release();

        program->setAttributeValue("colors",facet_color);
        vaos[0]->release();
        program->release();

    }
    //vao containing the data for the unselected lines
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[1]->bind();

        buffers[2].bind();
        buffers[2].allocate(positions_lines.data(), positions_lines.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[2].release();

        program->setAttributeValue("colors",edge_color);
        program->release();

        vaos[1]->release();

    }
    //vao containing the data for the points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[2]->bind();

        buffers[3].bind();
        buffers[3].allocate(positions_points.data(), positions_points.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[3].release();


        program->setAttributeValue("colors",vertex_color);
        buffers[6].release();
        program->release();

        vaos[2]->release();
    }

    are_buffers_filled = true;
}

void Scene_polyhedron_selection_item::compute_elements()
{
    positions_facets.clear();
    positions_lines.clear();
    positions_points.clear();
    normals.clear();
    //The facets
    {


        for(Selection_set_facet::iterator
            it = selected_facets.begin(),
            end = selected_facets.end();
            it != end; ++it)
        {
            const Kernel::Vector_3 n =
                    compute_facet_normal<Polyhedron::Facet,Kernel>(**it);

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());


            Polyhedron::Halfedge_around_facet_circulator
                    he = (*it)->facet_begin(),
                    cend = he;

            CGAL_For_all(he,cend)
            {
                const Kernel::Point_3& p = he->vertex()->point();
                positions_facets.push_back(p.x());
                positions_facets.push_back(p.y());
                positions_facets.push_back(p.z());
            }
        }
    }

    //The Lines
    {

        for(Selection_set_edge::iterator it = selected_edges.begin(); it != selected_edges.end(); ++it) {
            const Kernel::Point_3& a = (it->halfedge())->vertex()->point();
            const Kernel::Point_3& b = (it->halfedge())->opposite()->vertex()->point();
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
        }

    }

    //The points
    {
        for(Selection_set_vertex::iterator
            it = selected_vertices.begin(),
            end = selected_vertices.end();
            it != end; ++it)
        {
            const Kernel::Point_3& p = (*it)->point();
            positions_points.push_back(p.x());
            positions_points.push_back(p.y());
            positions_points.push_back(p.z());
        }
    }
}

void Scene_polyhedron_selection_item::draw(Viewer_interface* viewer) const
{

    if(!are_buffers_filled)
        initialize_buffers(viewer);

    draw_points(viewer);
    GLfloat offset_factor;
    GLfloat offset_units;
    qFunc.glGetFloatv( GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    qFunc.glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
    glPolygonOffset(-1.f, 1.f);

    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_TRIANGLES, 0, positions_facets.size()/3);
    program->release();
    vaos[0]->release();
    glPolygonOffset(offset_factor, offset_units);
    draw_edges(viewer);


}

void Scene_polyhedron_selection_item::draw_edges(Viewer_interface* viewer) const
{

    qFunc.glLineWidth(3.f);
    vaos[1]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
    program->release();
    vaos[1]->release();
    qFunc.glLineWidth(1.f);
}

void Scene_polyhedron_selection_item::draw_points(Viewer_interface* viewer) const
{
    qFunc.glPointSize(5.f);
    vaos[2]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_POINTS, 0, positions_points.size()/3);
    program->release();
    vaos[2]->release();
    qFunc.glPointSize(1.f);

}
