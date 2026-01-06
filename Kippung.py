import math
import tkinter as tk
from tkinter import ttk, messagebox


B_DEFAULT = 0.75
H_DEFAULT = 2.5
G_DEFAULT = 9.81


def calc_f_crit(m, h, B=B_DEFAULT, g=G_DEFAULT):
    a = B / 2.0
    return (m * g * a) / h


def calc_jet_force_from_Q(Q, rho, D0=None, A0=None, x=0.0, theta_deg=11.0, eta=1.0):
    if A0 is None:
        if D0 is None:
            raise ValueError("D0 or A0 must be provided.")
        A0 = (math.pi / 4.0) * (D0 ** 2)
    theta_rad = math.radians(theta_deg)
    D_x = math.sqrt(4.0 * A0 / math.pi) + 2.0 * x * math.tan(theta_rad)
    A_x = (math.pi / 4.0) * (D_x ** 2)
    v_x = Q / A_x
    # Very simplified model; real jet behavior, turbulence, and reflections can differ a lot.
    # Larger distance x reduces v_x due to increased jet area.
    return eta * rho * Q * v_x


def calc_force_from_thrust(T, eta=1.0):
    # Thrust is assumed to be the resulting force; distance is not applied again.
    return eta * T


def calc_wind_force(v_wind, A_wall, rho=1.225, Cd=1.2, gust_factor=1.0):
    q = 0.5 * rho * (v_wind ** 2)
    return q * Cd * A_wall * gust_factor


def _parse_float(value, name, default=None):
    raw = value.strip()
    if raw == "":
        if default is None:
            raise ValueError(f"{name} fehlt.")
        return default
    try:
        return float(raw)
    except ValueError as exc:
        raise ValueError(f"{name} ungueltig.") from exc


def _maybe_float(value):
    if value is None:
        return None
    raw = value.strip()
    if raw == "":
        return None
    try:
        return float(raw)
    except ValueError:
        raise ValueError("Ungueltige Eingabe.")


def _require_positive(val, name, allow_zero=False):
    if allow_zero and val == 0:
        return
    if val <= 0:
        raise ValueError(f"{name} muss > 0 sein.")


def _require_eta(val, name="eta"):
    if not (0 < val <= 1.0):
        raise ValueError(f"{name} muss in (0..1] liegen.")


class KippungApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Kippstabilitaet der Sichtschutzwand")
        self.root.geometry("980x620")
        self.root.minsize(920, 560)

        self._init_style()
        self._init_vars()
        self._build_layout()

    def _init_style(self):
        style = ttk.Style()
        style.theme_use("clam")
        base_font = ("Segoe UI", 10)
        title_font = ("Segoe UI Semibold", 14)

        self.root.configure(bg="#f5f1e8")
        style.configure("TFrame", background="#f5f1e8")
        style.configure("Card.TFrame", background="#ffffff", relief="flat")
        style.configure("TLabel", background="#f5f1e8", font=base_font)
        style.configure("Card.TLabel", background="#ffffff", font=base_font)
        style.configure("Title.TLabel", background="#f5f1e8", font=title_font, foreground="#2b2b2b")
        style.configure("TButton", font=base_font, padding=6)
        style.configure("Accent.TButton", font=base_font, padding=6, background="#1d6a5c", foreground="#ffffff")
        style.map("Accent.TButton", background=[("active", "#175548")])
        style.configure("TLabelframe", background="#f5f1e8", font=base_font)
        style.configure("TLabelframe.Label", background="#f5f1e8", font=("Segoe UI Semibold", 10))
        style.configure("TEntry", padding=4)

    def _init_vars(self):
        self.m_var = tk.StringVar()
        self.h_var = tk.StringVar()
        self.B_var = tk.StringVar(value=str(B_DEFAULT))
        self.H_var = tk.StringVar(value=str(H_DEFAULT))
        self.method_var = tk.StringVar(value="1")

        self.direct_F_var = tk.StringVar()

        self.jet_rho_var = tk.StringVar(value="1.225")
        self.jet_Q_var = tk.StringVar()
        self.jet_geom_var = tk.StringVar(value="a")
        self.jet_D0_var = tk.StringVar()
        self.jet_A0_var = tk.StringVar()
        self.jet_x_var = tk.StringVar(value="0")
        self.jet_theta_var = tk.StringVar(value="11")
        self.jet_eta_var = tk.StringVar(value="1.0")

        self.thrust_T_var = tk.StringVar()
        self.thrust_eta_var = tk.StringVar(value="1.0")

        self.wind_rho_var = tk.StringVar(value="1.225")
        self.wind_v_var = tk.StringVar()
        self.wind_A_var = tk.StringVar()
        self.wind_Cd_var = tk.StringVar(value="1.2")
        self.wind_gust_var = tk.StringVar(value="1.0")
        self.wind_use_mid_var = tk.BooleanVar(value=False)

        self.result_var = tk.StringVar(value="Bitte Eingaben machen und berechnen.")

    def _build_layout(self):
        header = ttk.Frame(self.root)
        header.pack(fill="x", padx=18, pady=(16, 8))
        ttk.Label(header, text="Seitliche Kippstabilitaet", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            header,
            text="Berechnung fuer mobile Sichtschutzwand auf Rollen (starrer Koerper, vereinfachtes Modell).",
        ).pack(anchor="w", pady=(4, 0))

        body = ttk.Frame(self.root)
        body.pack(fill="both", expand=True, padx=18, pady=(0, 16))

        left = ttk.Frame(body)
        left.pack(side="left", fill="both", expand=True, padx=(0, 12))

        right = ttk.Frame(body)
        right.pack(side="right", fill="both", expand=True)

        self._build_general_panel(left)
        self._build_method_panel(left)
        self._build_action_panel(left)
        self._build_result_panel(right)

    def _build_general_panel(self, parent):
        frame = ttk.Labelframe(parent, text="Basisdaten")
        frame.pack(fill="x", pady=(0, 12))

        grid = ttk.Frame(frame)
        grid.pack(fill="x", padx=12, pady=10)

        ttk.Label(grid, text="Masse m [kg]").grid(row=0, column=0, sticky="w")
        ttk.Entry(grid, textvariable=self.m_var, width=16).grid(row=0, column=1, padx=(10, 20))

        ttk.Label(grid, text="Angriffshoehe h [m]").grid(row=0, column=2, sticky="w")
        ttk.Entry(grid, textvariable=self.h_var, width=16).grid(row=0, column=3, padx=(10, 0))

        ttk.Label(grid, text="Wandhoehe H [m]").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(grid, textvariable=self.H_var, width=16).grid(row=1, column=1, padx=(10, 20), pady=(8, 0))

        ttk.Label(grid, text="Standbreite B [m]").grid(row=1, column=2, sticky="w", pady=(8, 0))
        ttk.Entry(grid, textvariable=self.B_var, width=16).grid(row=1, column=3, padx=(10, 0), pady=(8, 0))

        ttk.Label(grid, text=f"g = {G_DEFAULT} m/s^2").grid(row=2, column=3, sticky="w", pady=(8, 0))

    def _build_method_panel(self, parent):
        frame = ttk.Labelframe(parent, text="Kraftmodell")
        frame.pack(fill="both", expand=True, pady=(0, 12))

        selector = ttk.Frame(frame)
        selector.pack(fill="x", padx=12, pady=8)

        ttk.Radiobutton(selector, text="Direkte Kraft", variable=self.method_var, value="1",
                        command=self._update_method_view).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(selector, text="Impeller/Luftstrahl", variable=self.method_var, value="2",
                        command=self._update_method_view).grid(row=0, column=1, sticky="w", padx=(12, 0))
        ttk.Radiobutton(selector, text="Impeller-Schub", variable=self.method_var, value="3",
                        command=self._update_method_view).grid(row=0, column=2, sticky="w", padx=(12, 0))
        ttk.Radiobutton(selector, text="Windlast", variable=self.method_var, value="4",
                        command=self._update_method_view).grid(row=0, column=3, sticky="w", padx=(12, 0))

        self.method_container = ttk.Frame(frame)
        self.method_container.pack(fill="both", expand=True, padx=12, pady=(0, 12))

        self.direct_frame = self._build_direct_frame(self.method_container)
        self.jet_frame = self._build_jet_frame(self.method_container)
        self.thrust_frame = self._build_thrust_frame(self.method_container)
        self.wind_frame = self._build_wind_frame(self.method_container)

        self._update_method_view()

    def _build_direct_frame(self, parent):
        frame = ttk.Frame(parent)
        ttk.Label(frame, text="Direkte Kraft F_applied [N]").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.direct_F_var, width=16).grid(row=0, column=1, padx=(10, 0))
        return frame

    def _build_jet_frame(self, parent):
        frame = ttk.Frame(parent)
        ttk.Label(frame, text="Luftdichte rho [kg/m^3]").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.jet_rho_var, width=16).grid(row=0, column=1, padx=(10, 20))

        ttk.Label(frame, text="Volumenstrom Q [m^3/s]").grid(row=0, column=2, sticky="w")
        ttk.Entry(frame, textvariable=self.jet_Q_var, width=16).grid(row=0, column=3, padx=(10, 0))

        geom = ttk.LabelFrame(frame, text="Auslassgeometrie")
        geom.grid(row=1, column=0, columnspan=4, sticky="ew", pady=(10, 0))

        ttk.Radiobutton(geom, text="Durchmesser D0 [m]", variable=self.jet_geom_var, value="a",
                        command=self._update_jet_geom).grid(row=0, column=0, sticky="w", padx=6, pady=4)
        ttk.Entry(geom, textvariable=self.jet_D0_var, width=14).grid(row=0, column=1, padx=(4, 20))

        ttk.Radiobutton(geom, text="Flaeche A0 [m^2]", variable=self.jet_geom_var, value="b",
                        command=self._update_jet_geom).grid(row=0, column=2, sticky="w", padx=6, pady=4)
        ttk.Entry(geom, textvariable=self.jet_A0_var, width=14).grid(row=0, column=3, padx=(4, 0))

        ttk.Label(frame, text="Abstand x [m]").grid(row=2, column=0, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.jet_x_var, width=16).grid(row=2, column=1, padx=(10, 20), pady=(10, 0))

        ttk.Label(frame, text="Aufweitungswinkel theta [Grad]").grid(row=2, column=2, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.jet_theta_var, width=16).grid(row=2, column=3, padx=(10, 0), pady=(10, 0))

        ttk.Label(frame, text="Impulsfaktor eta (0..1]").grid(row=3, column=0, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.jet_eta_var, width=16).grid(row=3, column=1, padx=(10, 0), pady=(10, 0))

        note = (
            "Hinweis: Vereinfachtes Freistrahl-Modell. Reale Effekte und Abprall koennen stark abweichen."
        )
        ttk.Label(frame, text=note, wraplength=520).grid(row=4, column=0, columnspan=4, sticky="w", pady=(10, 0))

        return frame

    def _build_thrust_frame(self, parent):
        frame = ttk.Frame(parent)
        ttk.Label(frame, text="Schub T [N]").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.thrust_T_var, width=16).grid(row=0, column=1, padx=(10, 20))

        ttk.Label(frame, text="Wirkungsgrad eta (0..1]").grid(row=0, column=2, sticky="w")
        ttk.Entry(frame, textvariable=self.thrust_eta_var, width=16).grid(row=0, column=3, padx=(10, 0))

        note = "Hinweis: T gilt bereits als resultierende Kraft; Abstand wird nicht zusaetzlich beruecksichtigt."
        ttk.Label(frame, text=note, wraplength=520).grid(row=1, column=0, columnspan=4, sticky="w", pady=(10, 0))

        return frame

    def _build_wind_frame(self, parent):
        frame = ttk.Frame(parent)
        ttk.Label(frame, text="Luftdichte rho [kg/m^3]").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.wind_rho_var, width=16).grid(row=0, column=1, padx=(10, 20))

        ttk.Label(frame, text="Windgeschwindigkeit v_wind [m/s]").grid(row=0, column=2, sticky="w")
        ttk.Entry(frame, textvariable=self.wind_v_var, width=16).grid(row=0, column=3, padx=(10, 0))

        ttk.Label(frame, text="Wandflaeche A_wall [m^2]").grid(row=1, column=0, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.wind_A_var, width=16).grid(row=1, column=1, padx=(10, 20), pady=(10, 0))

        ttk.Label(frame, text="Widerstandsbeiwert Cd").grid(row=1, column=2, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.wind_Cd_var, width=16).grid(row=1, column=3, padx=(10, 0), pady=(10, 0))

        ttk.Label(frame, text="Boeenfaktor gust_factor").grid(row=2, column=0, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.wind_gust_var, width=16).grid(row=2, column=1, padx=(10, 0), pady=(10, 0))

        ttk.Checkbutton(
            frame,
            text="Resultierende bei H/2 ansetzen (gleichmaessiger Druck)",
            variable=self.wind_use_mid_var,
        ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(10, 0))

        return frame

    def _build_action_panel(self, parent):
        frame = ttk.Frame(parent)
        frame.pack(fill="x")
        ttk.Button(frame, text="Berechnen", style="Accent.TButton", command=self._calculate).pack(
            side="left"
        )
        ttk.Button(frame, text="Zuruecksetzen", command=self._reset).pack(side="left", padx=(10, 0))
        ttk.Button(frame, text="Beenden", command=self.root.destroy).pack(side="right")

    def _build_result_panel(self, parent):
        card = ttk.Frame(parent, style="Card.TFrame")
        card.pack(fill="both", expand=True)
        inner = ttk.Frame(card, style="Card.TFrame")
        inner.pack(fill="both", expand=True, padx=16, pady=16)

        ttk.Label(inner, text="Ergebnis", font=("Segoe UI Semibold", 11), style="Card.TLabel").pack(
            anchor="w", pady=(0, 8)
        )

        self.result_label = ttk.Label(
            inner,
            textvariable=self.result_var,
            justify="left",
            style="Card.TLabel",
            wraplength=360,
        )
        self.result_label.pack(anchor="nw", fill="both", expand=True)

    def _update_method_view(self):
        for child in self.method_container.winfo_children():
            child.pack_forget()
        method = self.method_var.get()
        if method == "1":
            self.direct_frame.pack(fill="x")
        elif method == "2":
            self.jet_frame.pack(fill="x")
        elif method == "3":
            self.thrust_frame.pack(fill="x")
        else:
            self.wind_frame.pack(fill="x")

    def _update_jet_geom(self):
        geom = self.jet_geom_var.get()
        if geom == "a":
            self.jet_A0_var.set("")
        else:
            self.jet_D0_var.set("")

    def _calculate(self):
        try:
            m = _parse_float(self.m_var.get(), "Masse m")
            _require_positive(m, "Masse m")

            h = _parse_float(self.h_var.get(), "Angriffshoehe h")
            _require_positive(h, "Angriffshoehe h")

            H = _parse_float(self.H_var.get(), "Wandhoehe H", default=H_DEFAULT)
            _require_positive(H, "Wandhoehe H")

            B = _parse_float(self.B_var.get(), "Standbreite B", default=B_DEFAULT)
            _require_positive(B, "Standbreite B")

            if h > H:
                messagebox.showwarning(
                    "Hinweis", "Angriffspunkt liegt ueber der Wandhoehe."
                )

            method = self.method_var.get()
            F_applied = None
            method_name = "nicht angegeben"
            info_note = None
            if method == "1":
                F_val = _maybe_float(self.direct_F_var.get())
                if F_val is None:
                    info_note = "Keine direkte Kraft eingetragen. Es wird nur F_crit berechnet."
                else:
                    _require_positive(F_val, "Direkte Kraft F")
                    F_applied = F_val
                    method_name = "Direkte Kraft"
                h_used = h
            elif method == "2":
                Q_val = _maybe_float(self.jet_Q_var.get())
                D0_val = _maybe_float(self.jet_D0_var.get())
                A0_val = _maybe_float(self.jet_A0_var.get())
                if Q_val is None or (self.jet_geom_var.get() == "a" and D0_val is None) or \
                        (self.jet_geom_var.get() == "b" and A0_val is None):
                    info_note = "Unvollstaendige Jet-Eingaben. Es wird nur F_crit berechnet."
                else:
                    rho = _parse_float(self.jet_rho_var.get(), "rho", default=1.225)
                    _require_positive(rho, "rho")
                    Q = _parse_float(self.jet_Q_var.get(), "Volumenstrom Q")
                    _require_positive(Q, "Volumenstrom Q")
                    geom = self.jet_geom_var.get()
                    D0 = None
                    A0 = None
                    if geom == "a":
                        D0 = _parse_float(self.jet_D0_var.get(), "Durchmesser D0")
                        _require_positive(D0, "Durchmesser D0")
                    else:
                        A0 = _parse_float(self.jet_A0_var.get(), "Flaeche A0")
                        _require_positive(A0, "Flaeche A0")
                    x = _parse_float(self.jet_x_var.get(), "Abstand x", default=0.0)
                    _require_positive(x, "Abstand x", allow_zero=True)
                    theta = _parse_float(self.jet_theta_var.get(), "theta", default=11.0)
                    _require_positive(theta, "theta")
                    eta = _parse_float(self.jet_eta_var.get(), "eta", default=1.0)
                    _require_eta(eta)
                    F_applied = calc_jet_force_from_Q(Q, rho, D0=D0, A0=A0, x=x, theta_deg=theta, eta=eta)
                    method_name = "Impeller/Luftstrahl"
                h_used = h
            elif method == "3":
                T_val = _maybe_float(self.thrust_T_var.get())
                if T_val is None:
                    info_note = "Kein Schub eingetragen. Es wird nur F_crit berechnet."
                else:
                    _require_positive(T_val, "Schub T")
                    eta = _parse_float(self.thrust_eta_var.get(), "eta", default=1.0)
                    _require_eta(eta)
                    F_applied = calc_force_from_thrust(T_val, eta=eta)
                    method_name = "Impeller-Schub"
                h_used = h
            else:
                v_val = _maybe_float(self.wind_v_var.get())
                A_val = _maybe_float(self.wind_A_var.get())
                if v_val is None or A_val is None:
                    info_note = "Unvollstaendige Wind-Eingaben. Es wird nur F_crit berechnet."
                else:
                    rho = _parse_float(self.wind_rho_var.get(), "rho", default=1.225)
                    _require_positive(rho, "rho")
                    v_wind = _parse_float(self.wind_v_var.get(), "Windgeschwindigkeit v_wind")
                    _require_positive(v_wind, "Windgeschwindigkeit v_wind")
                    A_wall = _parse_float(self.wind_A_var.get(), "Wandflaeche A_wall")
                    _require_positive(A_wall, "Wandflaeche A_wall")
                    Cd = _parse_float(self.wind_Cd_var.get(), "Cd", default=1.2)
                    _require_positive(Cd, "Cd")
                    gust = _parse_float(self.wind_gust_var.get(), "gust_factor", default=1.0)
                    _require_positive(gust, "gust_factor")
                    F_applied = calc_wind_force(v_wind, A_wall, rho=rho, Cd=Cd, gust_factor=gust)
                    method_name = "Windlast"
                if self.wind_use_mid_var.get():
                    h_used = H / 2.0
                else:
                    h_used = h

            self._update_result(m, h_used, F_applied, method_name, H, B, info_note=info_note)
        except ValueError as exc:
            messagebox.showerror("Eingabefehler", str(exc))

    def _update_result(self, m, h_used, F_applied, method_name, H, B, info_note=None):
        F_crit = calc_f_crit(m, h_used, B=B, g=G_DEFAULT)
        if F_applied is not None:
            sf = F_crit / F_applied if F_applied > 0 else float("inf")
        else:
            sf = None

        if F_applied is None:
            status = "Kein Vergleich moeglich, da keine Kraft angegeben wurde."
        else:
            status = "KIPPEN nach diesem Modell zu erwarten." if F_applied >= F_crit else \
                "standsicher bzgl. Kippen nach diesem Modell."

        text = f"F_crit = {F_crit:.3f} N\n"

        if F_applied is None:
            text += (
                "F_applied = - (nicht angegeben)\n"
                "Sicherheitsfaktor SF = -\n\n"
            )
        else:
            text += (
                f"F_applied = {F_applied:.3f} N ({method_name})\n"
                f"Sicherheitsfaktor SF = {sf:.3f}\n\n"
            )

        if info_note:
            text += f"Hinweis: {info_note}\n\n"

        text += f"Bewertung: {status}"
        self.result_var.set(text)

    def _reset(self):
        self.m_var.set("")
        self.h_var.set("")
        self.H_var.set(str(H_DEFAULT))
        self.B_var.set(str(B_DEFAULT))
        self.method_var.set("1")
        self.direct_F_var.set("")
        self.jet_Q_var.set("")
        self.jet_D0_var.set("")
        self.jet_A0_var.set("")
        self.jet_x_var.set("0")
        self.jet_theta_var.set("11")
        self.jet_eta_var.set("1.0")
        self.thrust_T_var.set("")
        self.thrust_eta_var.set("1.0")
        self.wind_v_var.set("")
        self.wind_A_var.set("")
        self.wind_Cd_var.set("1.2")
        self.wind_gust_var.set("1.0")
        self.wind_use_mid_var.set(False)
        self.result_var.set("Bitte Eingaben machen und berechnen.")
        self._update_method_view()


def main():
    root = tk.Tk()
    KippungApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()

# Beispielwerte (kommentiert):
# m=100 kg, h=1.25 m, Wind: v=15 m/s, A=7.5 m^2, Cd=1.2, gust=1.3
# Jet: Q=0.8 m^3/s, D0=0.30 m, x=2.0 m, theta=11 Grad, eta=0.8
