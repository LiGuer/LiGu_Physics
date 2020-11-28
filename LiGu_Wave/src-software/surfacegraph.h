#ifndef SURFACEGRAPH_H
#define SURFACEGRAPH_H
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtWidgets/QSlider>
#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/Q3DTheme>
#include <QtGui/QImage>
using namespace QtDataVisualization;

class SurfaceGraph : public QObject
{
    Q_OBJECT
public:
    explicit SurfaceGraph(Q3DSurface *surface);
    ~SurfaceGraph();
    void fillDataProxy(const int N, double* Map);
    void initSurface();
    void ViewSurface(const int N, double* Map);

public Q_SLOTS:
private:
    Q3DSurface *m_graph = new Q3DSurface();
    QSurfaceDataProxy *m_DataProxy = new QSurfaceDataProxy();
    QSurface3DSeries *m_DataSeries = new QSurface3DSeries(m_DataProxy);
};

#endif // SURFACEGRAPH_H
