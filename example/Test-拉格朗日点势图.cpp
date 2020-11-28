int main()
{
	const double hx = 5E8;
	Universe univer;
	Dynamics Dyna;
	Graphics* g = new Graphics;
	Plot * p = new Plot;
	p->init(g);
	p->setAxisRange(-hx, -hx, hx, hx);
	Star* newStar;
	/*---------------- 地球 ----------------*/
	newStar = new Star;
	{
		double r[] = { 0,0 * hx,0 * hx };
		double v[] = { 0,0,0 };
		double m = 5.965E24;
		for (int i = 0; i < 3; i++) {
			newStar->r[i] = r[i];	newStar->v[i] = v[i];
		}
		newStar->mass = m;
	}univer.addStar(newStar);
	/*---------------- 月球 ----------------*/
	newStar = new Star;
	{
		double r[] = { 3.844039E8,0 * hx,0 * hx };
		double v[] = { 0,1E3,0 };
		double m = 7.349E22;
		for (int i = 0; i < 3; i++) {
			newStar->r[i] = r[i];	newStar->v[i] = v[i];
		}
		newStar->mass = m;
	}univer.addStar(newStar);

	/*---------------- 势图 ----------------*/
	Mat<double> map;
	map.zero(1000, 1000);
	Star stars[100];
	for (int k = 0; k < univer.StarNum; k++) {
		stars[k] = *univer.Stars[k];
	}
	Particle* parC = Dyna.MassCentre(stars, univer.StarNum);

	for (int y = 0; y < 1000; y++) {
		for (int x = 0; x < 1000; x++) {
			double E = 0;
			Vector posi = Vector(3, 1);
			posi[0] = p->pix2coor(x, 0); posi[1] = p->pix2coor(y, 1); posi[2] = 0;
			/*---------------- 引力势 ----------------*/
			E += Dyna.GravitatePotential(stars, univer.StarNum,posi);
			/*---------------- 惯性离心势 ----------------*/
			double R = univer.Stars[1]->r[0];
			double w = sqrt(Dyna.G * (univer.Stars[0]->mass + univer.Stars[1]->mass) / (R * R * R));
			E += Dyna.CentrifugalPotential(w, parC->r, posi);
			map[y * map.cols + x] = E;
		}
	}
	/*---------------- 遮盖 ----------------*/
	double Emax = -9999999;
	for (int y = 0; y < 1000; y++) {
		for (int x = 0; x < 1000; x++) {
			Emax = Emax >= map[y * map.cols + x] ? Emax : map[y * map.cols + x];
		}
	}
	double mint = 1.2 * Emax;
	for (int y = 0; y < 1000; y++) {
		for (int x = 0; x < 1000; x++) {
			if (map[y * map.cols + x] < mint)map[y * map.cols + x] = mint;
		}
	}
	g->PaintColor = 0xffffff;
	p->contourface(&map, 100);
	g->PicWrite("D:/LIGU.ppm");
}
