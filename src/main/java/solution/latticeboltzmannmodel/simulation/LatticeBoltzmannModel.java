package solution.latticeboltzmannmodel.simulation;

import javafx.application.Application;
import javafx.stage.Stage;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;
import javafx.scene.image.WritableImage;
import javafx.scene.image.ImageView;
import javafx.animation.AnimationTimer;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ChoiceBox;
import javafx.scene.paint.Color;
import javafx.geometry.Insets;
import javafx.geometry.Pos;

public class LatticeBoltzmannModel extends Application {
    private static final int NX_2D = 100;
    private static final int NY_2D = 100;
    private static final int NX_3D = 50;
    private static final int NY_3D = 50;
    private static final int NZ_3D = 50;
    private enum LatticeType { D2Q9, D3Q15, D3Q19, D3Q27 }
    private LatticeType currentLattice = LatticeType.D2Q9;
    private enum CollisionOperator { SRT, MRT, TRT, ENTROPIC }
    private CollisionOperator currentCollision = CollisionOperator.SRT;
    private enum BoundaryCondition { BOUNCE_BACK, VELOCITY, PRESSURE, PERIODIC, INLET_OUTLET, OPEN }
    private BoundaryCondition currentBoundary = BoundaryCondition.BOUNCE_BACK;
    private int Q;
    private int[][] C;
    private double[] W;
    private int NX, NY, NZ;
    /// D2Q9
    private static final int[][] C_D2Q9 = {
            {0, 0}, {1, 0}, {-1, 0},
            {0, 1}, {0, -1}, {1, 1},
            {-1, -1}, {-1, 1}, {1, -1}
    };
    private static final double[] W_D2Q9 = {
            4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    /// D3Q15
    private static final int[][] C_D3Q15 = {
            {0, 0, 0}, {1, 0, 0}, {-1, 0, 0},
            {0, 1, 0}, {0, -1, 0}, {0, 0, 1},
            {0, 0, -1}, {1, 1, 1}, {-1, -1, -1},
            {1, 1, -1}, {-1, -1, 1}, {1, -1, 1},
            {-1, 1, -1}, {1, -1, -1}, {-1, 1, 1}
    };
    private static final double[] W_D3Q15 = {
            2.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
            1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0
    };
    /// D3Q19
    private static final int[][] C_D3Q19 = {
            {0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {1, 1, 0}, {-1, -1, 0},
            {1, -1, 0}, {-1, 1, 0}, {1, 0, 1}, {-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, -1}, {0, 1, -1}, {0, -1, 1}
    };
    private static final double[] W_D3Q19 = {
            1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };

    /// D3Q27
    private static final int[][] C_D3Q27 = {
            {0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {1, 1, 0}, {-1, -1, 0},
            {1, -1, 0}, {-1, 1, 0}, {1, 0, 1}, {-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, -1}, {0, 1, -1},
            {0, -1, 1}, {1, 1, 1}, {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {1, -1, -1}, {-1, 1, 1}
    };
    private static final double[] W_D3Q27 = {
            8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
            1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0
    };
    private static final double TAU = 0.6;
    private static final double OMEGA = 1.0 / TAU;
    private static final double U0 = 0.1;
    private static final double RHO0 = 1.0;
    private static final double[] S_MRT_D2Q9 = {
            0.0, 1.4, 1.2, 0.0, 1.0, 0.0, 1.0, 1.0 / TAU, 1.0 / TAU
    };
    private static final double[][] M_D2Q9 = {
            {1, 1, 1, 1, 1, 1, 1, 1, 1},
            {-4, -1, -1, -1, -1, 2, 2, 2, 2},
            {4, -2, -2, -2, -2, 1, 1, 1, 1},
            {0, 1, -1, 0, 0, 1, -1, -1, 1},
            {0, 0, 0, 1, -1, 1, -1, 1, -1},
            {0, 2, 2, -1, -1, 0, 0, 0, 0},
            {0, 0, 0, 1, 1, -2, -2, 2, 2},
            {0, 1, -1, 0, 0, -1, 1, 1, -1},
            {0, 0, 0, 1, -1, -1, 1, -1, 1}
    };
    private static final double TAU_S = 0.6;
    private static final double TAU_A = 0.5;
    private double[][][][][] f; // [x][y][z][q][current/next]
    private int current = 0;
    private int next = 1;
    private double[][][] rho;
    private double[][][] ux;
    private double[][][] uy;
    private double[][][] uz;
    private boolean[][][] obstacle;
    private WritableImage image;
    private ImageView imageView;
    private AnimationTimer timer;
    private boolean isRunning = false;
    @Override
    public void start(Stage primaryStage) {
        int canvasSize = 360;
        image = new WritableImage(canvasSize, canvasSize);
        imageView = new ImageView(image);
        imageView.setFitWidth(canvasSize);
        imageView.setFitHeight(canvasSize);
        imageView.setPreserveRatio(true);
        initialize();
        /// Control panel
        Button startButton = new Button("Start");
        Button stopButton = new Button("Stop");
        Button resetButton = new Button("Reset");
        ChoiceBox<LatticeType> latticeSelector = new ChoiceBox<>();
        ChoiceBox<CollisionOperator> collisionSelector = new ChoiceBox<>();
        ChoiceBox<BoundaryCondition> boundarySelector = new ChoiceBox<>();
        latticeSelector.getItems().addAll(LatticeType.D2Q9, LatticeType.D3Q15, LatticeType.D3Q19, LatticeType.D3Q27);
        collisionSelector.getItems().addAll(CollisionOperator.SRT, CollisionOperator.MRT, CollisionOperator.TRT, CollisionOperator.ENTROPIC);
        boundarySelector.getItems().addAll(BoundaryCondition.BOUNCE_BACK, BoundaryCondition.VELOCITY, BoundaryCondition.PRESSURE, BoundaryCondition.PERIODIC, BoundaryCondition.INLET_OUTLET, BoundaryCondition.OPEN);
        latticeSelector.setValue(LatticeType.D2Q9);
        collisionSelector.setValue(CollisionOperator.SRT);
        boundarySelector.setValue(BoundaryCondition.BOUNCE_BACK);
        /// Style controls
        String buttonStyle = "-fx-background-color: #333333; -fx-text-fill: white; -fx-font-size: 14; -fx-font-weight: bold; " +
                "-fx-background-radius: 20; -fx-padding: 10 20; -fx-min-height: 40; -fx-effect: dropshadow(gaussian, rgba(0,0,0,0.3), 5, 0, 0, 2);";
        String choiceBoxStyle = "-fx-background-color: #333333; -fx-text-fill: white; -fx-font-size: 14; -fx-mark-color: white; " +
                "-fx-background-radius: 10; -fx-padding: 10; -fx-min-height: 40; -fx-effect: dropshadow(gaussian, rgba(0,0,0,0.3), 5, 0, 0, 2);";

        startButton.setStyle(buttonStyle);
        stopButton.setStyle(buttonStyle);
        resetButton.setStyle(buttonStyle);
        latticeSelector.setStyle(choiceBoxStyle);
        collisionSelector.setStyle(choiceBoxStyle);
        boundarySelector.setStyle(choiceBoxStyle);
        /// Button widths
        startButton.setPrefWidth(200);
        stopButton.setPrefWidth(200);
        resetButton.setPrefWidth(200);
        latticeSelector.setPrefWidth(200);
        collisionSelector.setPrefWidth(200);
        boundarySelector.setPrefWidth(200);
        /// Button actions
        startButton.setOnAction(e -> startSimulation());
        stopButton.setOnAction(e -> stopSimulation());
        resetButton.setOnAction(e -> resetSimulation());
        latticeSelector.setOnAction(e -> {
            currentLattice = latticeSelector.getValue();
            resetSimulation();
        });
        collisionSelector.setOnAction(e -> {
            currentCollision = collisionSelector.getValue();
            resetSimulation();
        });
        boundarySelector.setOnAction(e -> {
            currentBoundary = boundarySelector.getValue();
            resetSimulation();
        });
        /// Control panel layout
        Label titleLabel = new Label("Lattice Boltzmann Model");
        titleLabel.setStyle("-fx-text-fill: white; -fx-font-size: 18; -fx-font-weight: bold; -fx-padding: 10;");
        VBox controlPanel = new VBox(10, titleLabel, startButton, stopButton, resetButton, latticeSelector, collisionSelector, boundarySelector);
        controlPanel.setPadding(new Insets(10));
        controlPanel.setStyle("-fx-background-color: black;");
        controlPanel.setAlignment(Pos.CENTER);
        /// Main layout
        Pane simulationPane = new Pane(imageView);
        simulationPane.setStyle("-fx-background-color: black;");
        simulationPane.setPrefHeight(canvasSize);
        VBox mainLayout = new VBox(simulationPane, controlPanel);
        mainLayout.setStyle("-fx-background-color: black;");
        /// Set up
        Scene scene = new Scene(mainLayout, 360, 640);
        scene.setFill(Color.BLACK);
        primaryStage.setTitle("Lattice Boltzmann Model Simulation");
        primaryStage.setScene(scene);
        primaryStage.setResizable(false);
        primaryStage.show();
        /// Animation timer
        timer = new AnimationTimer() {
            @Override
            public void handle(long now) {
                if (isRunning) {
                    step();
                    visualize();
                }
            }
        };
    }
    private void initialize() {
        switch (currentLattice) {
            case D2Q9:
                Q = 9;
                C = C_D2Q9;
                W = W_D2Q9;
                NX = NX_2D;
                NY = NY_2D;
                NZ = 1;
                break;
            case D3Q15:
                Q = 15;
                C = C_D3Q15;
                W = W_D3Q15;
                NX = NX_3D;
                NY = NY_3D;
                NZ = NZ_3D;
                break;
            case D3Q19:
                Q = 19;
                C = C_D3Q19;
                W = W_D3Q19;
                NX = NX_3D;
                NY = NY_3D;
                NZ = NZ_3D;
                break;
            case D3Q27:
                Q = 27;
                C = C_D3Q27;
                W = W_D3Q27;
                NX = NX_3D;
                NY = NY_3D;
                NZ = NZ_3D;
                break;
        }
        f = new double[NX][NY][NZ][Q][2];
        rho = new double[NX][NY][NZ];
        ux = new double[NX][NY][NZ];
        uy = new double[NX][NY][NZ];
        uz = new double[NX][NY][NZ];
        obstacle = new boolean[NX][NY][NZ];
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    rho[x][y][z] = 0.0;
                    ux[x][y][z] = 0.0;
                    uy[x][y][z] = 0.0;
                    uz[x][y][z] = 0.0;
                    obstacle[x][y][z] = false;
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][0] = 0.0;
                        f[x][y][z][q][1] = 0.0;
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    double u0 = 0.0;
                    if (currentBoundary == BoundaryCondition.VELOCITY || currentBoundary == BoundaryCondition.INLET_OUTLET) {
                        u0 = x == 0 ? U0 : 0.0;
                    } else {
                        u0 = x < NX / 2 ? U0 : 0.0;
                    }
                    double[] feq = equilibrium(RHO0, u0, 0.0, 0.0);
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][current] = feq[q];
                        f[x][y][z][q][next] = feq[q];
                    }
                }
            }
        }
        if (currentLattice == LatticeType.D2Q9) {
            int cx = NX / 4;
            int cy = NY / 2;
            int radius = 10;
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    if (Math.sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy)) < radius) {
                        obstacle[x][y][0] = true;
                    }
                }
            }
        } else {
            int cx = NX / 4;
            int cy = NY / 2;
            int cz = NZ / 2;
            int radius = 5;
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int z = 0; z < NZ; z++) {
                        if (Math.sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz)) < radius) {
                            obstacle[x][y][z] = true;
                        }
                    }
                }
            }
        }
        current = 0;
        next = 1;
        computeMacroscopic();
        clearImage();
        visualize();
    }
    private void clearImage() {
        int canvasSize = 360;
        for (int x = 0; x < canvasSize; x++) {
            for (int y = 0; y < canvasSize; y++) {
                image.getPixelWriter().setArgb(x, y, 0xFF000000);
            }
        }
    }
    private void startSimulation() {
        if (!isRunning) {
            isRunning = true;
            timer.start();
        }
    }
    private void stopSimulation() {
        isRunning = false;
    }
    private void resetSimulation() {
        isRunning = false;
        timer.stop();
        initialize();
    }
    private double[] equilibrium(double rho, double ux, double uy, double uz) {
        double[] feq = new double[Q];
        double u2 = ux * ux + uy * uy + uz * uz;
        for (int q = 0; q < Q; q++) {
            double cu = C[q][0] * ux + C[q][1] * uy + (C[q].length > 2 ? C[q][2] * uz : 0.0);
            feq[q] = W[q] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
        }
        return feq;
    }
    private void computeMacroscopic() {
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    if (obstacle[x][y][z]) {
                        rho[x][y][z] = 0.0;
                        ux[x][y][z] = 0.0;
                        uy[x][y][z] = 0.0;
                        uz[x][y][z] = 0.0;
                        continue;
                    }
                    double density = 0.0;
                    double vx = 0.0;
                    double vy = 0.0;
                    double vz = 0.0;
                    for (int q = 0; q < Q; q++) {
                        density += f[x][y][z][q][current];
                        vx += f[x][y][z][q][current] * C[q][0];
                        vy += f[x][y][z][q][current] * C[q][1];
                        if (C[q].length > 2) vz += f[x][y][z][q][current] * C[q][2];
                    }
                    rho[x][y][z] = density;
                    ux[x][y][z] = density > 0 ? vx / density : 0.0;
                    uy[x][y][z] = density > 0 ? vy / density : 0.0;
                    uz[x][y][z] = density > 0 ? vz / density : 0.0;
                }
            }
        }
    }
    private void collide() {
        if (currentLattice != LatticeType.D2Q9 && currentCollision != CollisionOperator.SRT) {
            collideSRT();
            return;
        }
        switch (currentCollision) {
            case SRT:
                collideSRT();
                break;
            case MRT:
                collideMRT();
                break;
            case TRT:
                collideTRT();
                break;
            case ENTROPIC:
                collideEntropic();
                break;
        }
    }
    private void collideSRT() {
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    if (obstacle[x][y][z]) continue;
                    double[] feq = equilibrium(rho[x][y][z], ux[x][y][z], uy[x][y][z], uz[x][y][z]);
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][next] = f[x][y][z][q][current] + OMEGA * (feq[q] - f[x][y][z][q][current]);
                    }
                }
            }
        }
    }
    private void collideMRT() {
        double[][] invM = new double[9][9];
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 9; j++) {
                invM[i][j] = M_D2Q9[j][i] / (i == 0 ? 9.0 : (i == 1 || i == 2 ? 36.0 : 6.0));
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    if (obstacle[x][y][z]) continue;
                    double[] fi = new double[Q];
                    for (int q = 0; q < Q; q++) {
                        fi[q] = f[x][y][z][q][current];
                    }
                    double[] m = new double[Q];
                    for (int i = 0; i < Q; i++) {
                        m[i] = 0.0;
                        for (int j = 0; j < Q; j++) {
                            m[i] += M_D2Q9[i][j] * fi[j];
                        }
                    }
                    double[] meq = new double[Q];
                    meq[0] = rho[x][y][z];
                    meq[1] = -2.0 * rho[x][y][z] + 3.0 * rho[x][y][z] * (ux[x][y][z] * ux[x][y][z] + uy[x][y][z] * uy[x][y][z]);
                    meq[2] = rho[x][y][z] - 3.0 * rho[x][y][z] * (ux[x][y][z] * ux[x][y][z] + uy[x][y][z] * uy[x][y][z]);
                    meq[3] = rho[x][y][z] * ux[x][y][z];
                    meq[4] = rho[x][y][z] * uy[x][y][z];
                    meq[5] = rho[x][y][z] * (ux[x][y][z] * ux[x][y][z] - uy[x][y][z] * uy[x][y][z]);
                    meq[6] = rho[x][y][z] * ux[x][y][z] * uy[x][y][z];
                    meq[7] = 0.0;
                    meq[8] = 0.0;
                    for (int i = 0; i < Q; i++) {
                        m[i] -= S_MRT_D2Q9[i] * (m[i] - meq[i]);
                    }
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][next] = 0.0;
                        for (int i = 0; i < Q; i++) {
                            f[x][y][z][q][next] += invM[q][i] * m[i];
                        }
                    }
                }
            }
        }
    }
    private void collideTRT() {
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    if (obstacle[x][y][z]) continue;
                    double[] feq = equilibrium(rho[x][y][z], ux[x][y][z], uy[x][y][z], uz[x][y][z]);
                    for (int q = 0; q < Q; q++) {
                        int qOpp = getOppositeDirection(q);
                        double fs = 0.5 * (f[x][y][z][q][current] + f[x][y][z][qOpp][current]);
                        double fa = 0.5 * (f[x][y][z][q][current] - f[x][y][z][qOpp][current]);
                        double feqs = 0.5 * (feq[q] + feq[qOpp]);
                        double feqa = 0.5 * (feq[q] - feq[qOpp]);
                        f[x][y][z][q][next] = f[x][y][z][q][current] + (1.0 / TAU_S) * (feqs - fs) + (1.0 / TAU_A) * (feqa - fa);
                    }
                }
            }
        }
    }
    private void collideEntropic() {
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    if (obstacle[x][y][z]) continue;
                    double[] feq = equilibrium(rho[x][y][z], ux[x][y][z], uy[x][y][z], uz[x][y][z]);
                    double[] fi = new double[Q];
                    for (int q = 0; q < Q; q++) {
                        fi[q] = f[x][y][z][q][current];
                    }
                    double alpha = 1.0;
                    double[] fnew = new double[Q];
                    for (int iter = 0; iter < 10; iter++) {
                        double entropy = 0.0;
                        for (int q = 0; q < Q; q++) {
                            fnew[q] = fi[q] + alpha * (feq[q] - fi[q]);
                            if (fnew[q] > 0) entropy -= fnew[q] * Math.log(fnew[q] / W[q]);
                        }
                        double dEntropy = 0.0;
                        for (int q = 0; q < Q; q++) {
                            if (fnew[q] > 0) dEntropy -= (feq[q] - fi[q]) * (1.0 + Math.log(fnew[q] / W[q]));
                        }
                        if (Math.abs(dEntropy) < 1e-6) break;
                        alpha -= entropy / dEntropy;
                        if (alpha < 0) alpha = 0.5;
                    }
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][next] = fi[q] + alpha * (feq[q] - fi[q]);
                    }
                }
            }
        }
    }
    private void stream() {
        double[][][][] temp = new double[NX][NY][NZ][Q];
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    for (int q = 0; q < Q; q++) {
                        int xNext = x + C[q][0];
                        int yNext = y + C[q][1];
                        int zNext = (C[q].length > 2) ? z + C[q][2] : z;
                        if (currentBoundary == BoundaryCondition.PERIODIC) {
                            xNext = (xNext + NX) % NX;
                            yNext = (yNext + NY) % NY;
                            zNext = (zNext + NZ) % NZ;
                        } else {
                            if (xNext < 0 || xNext >= NX || yNext < 0 || yNext >= NY || zNext < 0 || zNext >= NZ) {
                                continue;
                            }
                        }
                        if (obstacle[xNext][yNext][zNext]) {
                            int qOpp = getOppositeDirection(q);
                            temp[x][y][z][qOpp] = f[x][y][z][q][next];
                        } else {
                            temp[xNext][yNext][zNext][q] = f[x][y][z][q][next];
                        }
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    for (int q = 0; q < Q; q++) {
                        f[x][y][z][q][current] = temp[x][y][z][q];
                    }
                }
            }
        }
    }
    private void applyBoundaryConditions() {
        switch (currentBoundary) {
            case BOUNCE_BACK:
                applyBounceBack();
                break;
            case VELOCITY:
                applyVelocity();
                break;
            case PRESSURE:
                applyPressure();
                break;
            case PERIODIC:
                break;
            case INLET_OUTLET:
                applyInletOutlet();
                break;
            case OPEN:
                applyOpen();
                break;
        }
    }
    private void applyBounceBack() {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][0] > 0) {
                        f[0][y][z][qOpp][current] = f[0][y][z][q][current];
                    }
                    if (C[q][0] < 0) {
                        f[NX-1][y][z][qOpp][current] = f[NX-1][y][z][q][current];
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][1] > 0) {
                        f[x][0][z][qOpp][current] = f[x][0][z][q][current];
                    }
                    if (C[q][1] < 0) {
                        f[x][NY-1][z][qOpp][current] = f[x][NY-1][z][q][current];
                    }
                }
            }
        }
        if (NZ > 1) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int q = 0; q < Q; q++) {
                        int qOpp = getOppositeDirection(q);
                        if (C[q].length > 2 && C[q][2] > 0) {
                            f[x][y][0][qOpp][current] = f[x][y][0][q][current];
                        }
                        if (C[q].length > 2 && C[q][2] < 0) {
                            f[x][y][NZ-1][qOpp][current] = f[x][y][NZ-1][q][current];
                        }
                    }
                }
            }
        }
    }
    private void applyVelocity() {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                if (obstacle[0][y][z]) continue;
                double rhoWall = 0.0;
                for (int q = 0; q < Q; q++) {
                    rhoWall += f[0][y][z][q][current];
                }
                double uxWall = U0;
                double uyWall = 0.0;
                double uzWall = 0.0;
                if (currentLattice == LatticeType.D2Q9) {
                    f[0][y][z][1][current] = f[0][y][z][2][current] + (2.0 * rhoWall * uxWall) / 3.0;
                    f[0][y][z][5][current] = f[0][y][z][6][current] + (rhoWall * uxWall) / 6.0 + (rhoWall * uyWall) / 6.0;
                    f[0][y][z][8][current] = f[0][y][z][7][current] + (rhoWall * uxWall) / 6.0 - (rhoWall * uyWall) / 6.0;
                } else {
                    double[] feq = equilibrium(rhoWall, uxWall, uyWall, uzWall);
                    for (int q = 0; q < Q; q++) {
                        if (C[q][0] > 0) {
                            f[0][y][z][q][current] = feq[q];
                        }
                    }
                }
            }
        }
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][0] < 0) {
                        f[NX-1][y][z][qOpp][current] = f[NX-1][y][z][q][current];
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][1] > 0) {
                        f[x][0][z][qOpp][current] = f[x][0][z][q][current];
                    }
                    if (C[q][1] < 0) {
                        f[x][NY-1][z][qOpp][current] = f[x][NY-1][z][q][current];
                    }
                }
            }
        }
        if (NZ > 1) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int q = 0; q < Q; q++) {
                        int qOpp = getOppositeDirection(q);
                        if (C[q].length > 2 && C[q][2] > 0) {
                            f[x][y][0][qOpp][current] = f[x][y][0][q][current];
                        }
                        if (C[q].length > 2 && C[q][2] < 0) {
                            f[x][y][NZ-1][qOpp][current] = f[x][y][NZ-1][q][current];
                        }
                    }
                }
            }
        }
    }
    private void applyPressure() {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                if (obstacle[NX-1][y][z]) continue;
                double rhoWall = RHO0;
                double uxWall = 0.0;
                double uyWall = 0.0;
                double uzWall = 0.0;
                if (currentLattice == LatticeType.D2Q9) {
                    uxWall = -1.0 + (f[NX-1][y][z][0][current] + f[NX-1][y][z][1][current] + f[NX-1][y][z][3][current] +
                            f[NX-1][y][z][4][current] + f[NX-1][y][z][5][current] + f[NX-1][y][z][8][current] +
                            2.0 * (f[NX-1][y][z][2][current] + f[NX-1][y][z][6][current] + f[NX-1][y][z][7][current])) / rhoWall;
                    f[NX-1][y][z][2][current] = f[NX-1][y][z][1][current] + (2.0 * rhoWall * uxWall) / 3.0;
                    f[NX-1][y][z][6][current] = f[NX-1][y][z][5][current] + (rhoWall * uxWall) / 6.0 - (rhoWall * uyWall) / 6.0;
                    f[NX-1][y][z][7][current] = f[NX-1][y][z][8][current] + (rhoWall * uxWall) / 6.0 + (rhoWall * uyWall) / 6.0;
                } else {
                    double[] feq = equilibrium(rhoWall, uxWall, uyWall, uzWall);
                    for (int q = 0; q < Q; q++) {
                        if (C[q][0] < 0) {
                            f[NX-1][y][z][q][current] = feq[q];
                        }
                    }
                }
            }
        }
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][0] > 0) {
                        f[0][y][z][qOpp][current] = f[0][y][z][q][current];
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][1] > 0) {
                        f[x][0][z][qOpp][current] = f[x][0][z][q][current];
                    }
                    if (C[q][1] < 0) {
                        f[x][NY-1][z][qOpp][current] = f[x][NY-1][z][q][current];
                    }
                }
            }
        }
        if (NZ > 1) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int q = 0; q < Q; q++) {
                        int qOpp = getOppositeDirection(q);
                        if (C[q].length > 2 && C[q][2] > 0) {
                            f[x][y][0][qOpp][current] = f[x][y][0][q][current];
                        }
                        if (C[q].length > 2 && C[q][2] < 0) {
                            f[x][y][NZ-1][qOpp][current] = f[x][y][NZ-1][q][current];
                        }
                    }
                }
            }
        }
    }
    private void applyInletOutlet() {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                if (obstacle[0][y][z]) continue;
                double rhoWall = 0.0;
                for (int q = 0; q < Q; q++) {
                    rhoWall += f[0][y][z][q][current];
                }
                double uxWall = U0;
                double uyWall = 0.0;
                double uzWall = 0.0;
                if (currentLattice == LatticeType.D2Q9) {
                    f[0][y][z][1][current] = f[0][y][z][2][current] + (2.0 * rhoWall * uxWall) / 3.0;
                    f[0][y][z][5][current] = f[0][y][z][6][current] + (rhoWall * uxWall) / 6.0 + (rhoWall * uyWall) / 6.0;
                    f[0][y][z][8][current] = f[0][y][z][7][current] + (rhoWall * uxWall) / 6.0 - (rhoWall * uyWall) / 6.0;
                } else {
                    double[] feq = equilibrium(rhoWall, uxWall, uyWall, uzWall);
                    for (int q = 0; q < Q; q++) {
                        if (C[q][0] > 0) {
                            f[0][y][z][q][current] = feq[q];
                        }
                    }
                }
            }
        }
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                if (obstacle[NX-1][y][z]) continue;
                double rhoWall = RHO0;
                double uxWall = 0.0;
                double uyWall = 0.0;
                double uzWall = 0.0;
                if (currentLattice == LatticeType.D2Q9) {
                    uxWall = -1.0 + (f[NX-1][y][z][0][current] + f[NX-1][y][z][1][current] + f[NX-1][y][z][3][current] +
                            f[NX-1][y][z][4][current] + f[NX-1][y][z][5][current] + f[NX-1][y][z][8][current] +
                            2.0 * (f[NX-1][y][z][2][current] + f[NX-1][y][z][6][current] + f[NX-1][y][z][7][current])) / rhoWall;
                    f[NX-1][y][z][2][current] = f[NX-1][y][z][1][current] + (2.0 * rhoWall * uxWall) / 3.0;
                    f[NX-1][y][z][6][current] = f[NX-1][y][z][5][current] + (rhoWall * uxWall) / 6.0 - (rhoWall * uyWall) / 6.0;
                    f[NX-1][y][z][7][current] = f[NX-1][y][z][8][current] + (rhoWall * uxWall) / 6.0 + (rhoWall * uyWall) / 6.0;
                } else {
                    double[] feq = equilibrium(rhoWall, uxWall, uyWall, uzWall);
                    for (int q = 0; q < Q; q++) {
                        if (C[q][0] < 0) {
                            f[NX-1][y][z][q][current] = feq[q];
                        }
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    int qOpp = getOppositeDirection(q);
                    if (C[q][1] > 0) {
                        f[x][0][z][qOpp][current] = f[x][0][z][q][current];
                    }
                    if (C[q][1] < 0) {
                        f[x][NY-1][z][qOpp][current] = f[x][NY-1][z][q][current];
                    }
                }
            }
        }
        if (NZ > 1) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int q = 0; q < Q; q++) {
                        int qOpp = getOppositeDirection(q);
                        if (C[q].length > 2 && C[q][2] > 0) {
                            f[x][y][0][qOpp][current] = f[x][y][0][q][current];
                        }
                        if (C[q].length > 2 && C[q][2] < 0) {
                            f[x][y][NZ-1][qOpp][current] = f[x][y][NZ-1][q][current];
                        }
                    }
                }
            }
        }
    }
    private void applyOpen() {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    if (C[q][0] > 0 && !obstacle[0][y][z]) {
                        f[0][y][z][q][current] = f[1][y][z][q][current];
                    }
                    if (C[q][0] < 0 && !obstacle[NX-1][y][z]) {
                        f[NX-1][y][z][q][current] = f[NX-2][y][z][q][current];
                    }
                }
            }
        }
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int q = 0; q < Q; q++) {
                    if (C[q][1] > 0 && !obstacle[x][0][z]) {
                        f[x][0][z][q][current] = f[x][1][z][q][current];
                    }
                    if (C[q][1] < 0 && !obstacle[x][NY-1][z]) {
                        f[x][NY-1][z][q][current] = f[x][NY-2][z][q][current];
                    }
                }
            }
        }
        if (NZ > 1) {
            for (int x = 0; x < NX; x++) {
                for (int y = 0; y < NY; y++) {
                    for (int q = 0; q < Q; q++) {
                        if (C[q].length > 2 && C[q][2] > 0 && !obstacle[x][y][0]) {
                            f[x][y][0][q][current] = f[x][y][1][q][current];
                        }
                        if (C[q].length > 2 && C[q][2] < 0 && !obstacle[x][y][NZ-1]) {
                            f[x][y][NZ-1][q][current] = f[x][y][NZ-2][q][current];
                        }
                    }
                }
            }
        }
    }
    private int getOppositeDirection(int q) {
        if (currentLattice == LatticeType.D2Q9) {
            switch (q) {
                case 1: return 2;
                case 2: return 1;
                case 3: return 4;
                case 4: return 3;
                case 5: return 6;
                case 6: return 5;
                case 7: return 8;
                case 8: return 7;
                default: return 0;
            }
        } else {
            for (int i = 0; i < Q; i++) {
                if (C[i][0] == -C[q][0] && C[i][1] == -C[q][1] && (C[i].length < 3 || C[i][2] == -C[q][2])) {
                    return i;
                }
            }
            return 0;
        }
    }
    private void step() {
        computeMacroscopic();
        collide();
        stream();
        applyBoundaryConditions();
    }
    private void visualize() {
        int visNX = currentLattice == LatticeType.D2Q9 ? NX_2D : NX_3D;
        int visNY = currentLattice == LatticeType.D2Q9 ? NY_2D : NY_3D;
        int canvasSize = 360;
        double scale = (double) canvasSize / visNX;
        for (int x = 0; x < visNX; x++) {
            for (int y = 0; y < visNY; y++) {
                int z = 0;
                int pixelXStart = (int) (x * scale);
                int pixelYStart = (int) (y * scale);
                int pixelXEnd = (int) ((x + 1) * scale);
                int pixelYEnd = (int) ((y + 1) * scale);
                int color;
                if (obstacle[x][y][z]) {
                    color = 0xFF000000;
                } else {
                    double vel = Math.sqrt(ux[x][y][z] * ux[x][y][z] + uy[x][y][z] * uy[x][y][z] + uz[x][y][z] * uz[x][y][z]);
                    double maxVel = 0.2;
                    double ratio = Math.min(vel / maxVel, 1.0);
                    int r = (int) (255 * ratio);
                    int b = (int) (255 * (1.0 - ratio));
                    color = (0xFF << 24) | (r << 16) | b;
                }
                for (int px = pixelXStart; px < pixelXEnd && px < canvasSize; px++) {
                    for (int py = pixelYStart; py < pixelYEnd && py < canvasSize; py++) {
                        image.getPixelWriter().setArgb(px, py, color);
                    }
                }
            }
        }
    }
    public static void main(String[] args) {
        launch(args);
    }
}