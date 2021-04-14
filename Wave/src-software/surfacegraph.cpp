#include "surfacegraph.h"
using namespace QtDataVisualization;
SurfaceGraph::SurfaceGraph(Q3DSurface *surface): m_graph(surface)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);
}
SurfaceGraph::~SurfaceGraph()
{
    delete m_graph;
}

void SurfaceGraph::fillDataProxy(const int N, double* Map)
{
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(N);
    for (int x = 0 ; x < N ; x++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(N);
        for (int y = 0; y < N; y++) {
            double z = *(Map + x*N +y);
            (*newRow)[y].setPosition(QVector3D(y, z, x));
        }
        *dataArray << newRow;
    }
    m_DataProxy->resetArray(dataArray);
}
void SurfaceGraph::initSurface()
{
    m_DataSeries->setDrawMode(QSurface3DSeries::DrawSurface);
    m_DataSeries->setFlatShadingEnabled(false);
    m_graph->axisX()->setLabelAutoRotation(30);
    m_graph->axisY()->setLabelAutoRotation(90);
    m_graph->axisZ()->setLabelAutoRotation(30);

    QLinearGradient gr;
    gr.setColorAt(0.5, QColor (20, 210 , 255, 255));
    gr.setColorAt(-0.5, QColor (20, 210 , 255, 255));
    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
}
void SurfaceGraph::ViewSurface(const int N, double* Map)
{
    fillDataProxy(N,Map);
    m_graph->axisX()->setRange(0, N);
    m_graph->axisY()->setRange(-1.0f, 1.0f);
    m_graph->axisZ()->setRange(0, N);
    m_graph->removeSeries(m_DataSeries);
    m_graph->addSeries(m_DataSeries);
}
