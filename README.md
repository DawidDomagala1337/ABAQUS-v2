# ABAQUS-v2

Skrypt `radon_exhalation.py` implementuje numeryczny model dyfuzji radonu na
bazie równań (1–10) z załączonych notatek (Nazaroff i Nero, 1988). Wystarczy
podać wymagane parametry fizyczne materiału, a program obliczy strumień
wywiewu radonu z sześciu powierzchni elementu i zwróci średni strumień
powierzchniowy \(J_b\).

## Jak równania zostały wprowadzone do kodu?

- **Równanie 1 (dyfuzja z rozpadem):** w solverze `solve_concentration` pole
  źródłowe obliczane jest jako \(\lambda Q \rho / \epsilon\), a przyrost
  \(C\) w każdym węźle wewnętrznym jest aktualizowany ze wzoru

  \[
  C_{i,j,k} = \frac{D\,\nabla^2 C + \lambda Q \rho / \epsilon}{\lambda + 2D\,(1/\Delta x^2 + 1/\Delta y^2 + 1/\Delta z^2)}.
  \]

- **Równania brzegowe 2a–2c (C = 0 na powierzchniach):** wartości brzegowe
  siatki są cały czas równe zero, co odpowiada przyjętej koncentracji w
  porach na powierzchni próbki.

- **Równania 8a–8c (strumień Ficka):** funkcja `face_fluxes` liczy
  \(J = -\epsilon D\,\partial C/\partial n\) metodą różnic skończonych na
  każdej z sześciu ścian.

- **Równanie 9 (średni strumień \(J_b\)):** funkcja `average_exhalation`
  uśrednia strumień po powierzchni całego elementu z użyciem pól poszczególnych
  ścian.

## Użycie

```bash
python radon_exhalation.py a b h D lambda_ Q rho epsilon [--nodes NX NY NZ] [--tol TOL] [--max-iter N] [--example]
```

### Gdzie wpisać parametry?

Parametry podstawiasz **bezpośrednio w miejscu nazw** `a b h D lambda_ Q rho epsilon` w linii wywołania:

1. Otwórz terminal w katalogu projektu (tam gdzie leży `radon_exhalation.py`).
2. Wpisz polecenie `python radon_exhalation.py` i **podaj po kolei liczby** dla:
   - `a`, `b`, `h` – połówki wymiarów elementu (grubość, szerokość, wysokość) w metrach.
   - `D` – współczynnik dyfuzji radonu w materiale [m²/s].
   - `lambda_` – stała rozpadu radonu [1/s].
   - `Q` – zawartość radu (czynnik emanacji) [Bq/kg].
   - `rho` – gęstość objętościowa materiału [kg/m³].
   - `epsilon` – porowatość (udział objętościowy przestrzeni porowej).
3. (Opcjonalnie) dopisz `--nodes NX NY NZ`, `--tol TOL` lub `--max-iter N`, jeżeli chcesz zmienić gęstość siatki lub kryterium zbieżności.

Możesz również użyć przełącznika `--example`, aby od razu wczytać zestaw
parametrów z publikowanego przykładu (a=0,3 m, b=0,15 m, h=0,3 m, D=8·10⁻⁷ m²/s,
λ=2,1·10⁻⁶ 1/s, Q=20 Bq/kg, ρ=2400 kg/m³, ε=0,25, siatka 5×5×5):

```
python radon_exhalation.py --example
```

W trybie `--example` możesz nadpisać pojedyncze wartości (np. tylko `--nodes`),
podając je po nazwie przełącznika.

Możesz podejrzeć opis parametrów wywołując pomoc:

```bash
python radon_exhalation.py --help
```

Gdzie:

- `a`, `b`, `h` – połowy wymiarów elementu (grubość, szerokość, wysokość) w metrach.
- `D` – współczynnik dyfuzji radonu w materiale [m²/s].
- `lambda_` – stała rozpadu radonu [1/s].
- `Q` – zawartość radu (czynnik emanacji) [Bq/kg].
- `rho` – gęstość objętościowa materiału [kg/m³].
- `epsilon` – porowatość (udział objętościowy przestrzeni porowej) z zakresu 0,005–1.
- `--nodes` – liczba węzłów siatki w kierunkach \(x, y, z\); domyślnie 5 × 5 × 5.
- `--tol` – tolerancja zbieżności metody Gaussa–Seidla.
- `--max-iter` – maksymalna liczba iteracji solvera.

### Jak działa parametr `--nodes`?

- `--nodes NX NY NZ` ustawia **liczbę punktów siatki** w kierunkach \(x, y, z\).
- Minimalna wartość to **3 w każdym kierunku** – dwa węzły brzegowe plus co najmniej
  jeden węzeł wewnętrzny (kod wymusza ten warunek i zgłosi błąd poniżej 3).
- Większa liczba węzłów zwiększa dokładność rozwiązania, ale też czas i pamięć;
  złożoność rośnie wprost proporcjonalnie do \(N_x N_y N_z\).
- Praktycznie: 20–60 węzłów na wymiar daje dobre przybliżenia na laptopie;
  siatki rzędu 100³ mogą być już bardzo wolne.
- Rozmiar kroku siatki wynika z wymiarów geometrycznych: \(\Delta x = 2a/(N_x-1)\),
  analogicznie dla \(\Delta y\) i \(\Delta z\).

Przykład (publikowany zestaw danych):

```bash
python radon_exhalation.py --example
```

Wynik zawiera długość dyfuzji \(l = \sqrt{D / \lambda}\), średnie strumienie
na każdej powierzchni oraz średni strumień wywiewu \(J_b\) uśredniony po
powierzchni całego elementu. Dla przykładowego zestawu danych:

```
Radon exhalation model (Eqs. 1–10)
Diffusion length l = 6.1721e-01 m
Mean outward fluxes on faces [Bq m^-2 s^-1]:
  x-: 5.0614e-03
  x+: 5.0614e-03
  y-: 7.0351e-03
  y+: 7.0351e-03
  z-: 4.8361e-03
  z+: 4.8361e-03
Area-weighted mean exhalation J_b = 5.9919e-03 Bq m^-2 s^-1
```

## Interfejs webowy (NiceGUI)

Interfejs przeglądarkowy znajduje się w oddzielnych plikach `nicegui_app.py`
oraz `FrontAPP.py` i korzysta z tych samych funkcji solvera co wersja CLI.
`nicegui_app.py` zawiera pełny formularz wraz z przyciskiem ładowania
parametrów przykładowych, natomiast `FrontAPP.py` jest lekkim oknem wejściowym
uruchamianym na porcie `8001` z podglądem postępu obliczeń. Aby skorzystać z
aplikacji:

1. Zainstaluj NiceGUI (jednorazowo):

   ```bash
   pip install nicegui
   ```

2. Wystartuj jedną z aplikacji:

   ```bash
   # pełny formularz z ładowaniem przykładu (port 8080)
   python nicegui_app.py

   # uproszczone okno wejściowe (port 8001)
   python FrontAPP.py
   ```

3. Otwórz w przeglądarce <http://localhost:8080> dla `nicegui_app.py` lub
   <http://127.0.0.1:8001> dla `FrontAPP.py`. W obu przypadkach wprowadź
   parametry modelu i kliknij przycisk obliczeń, aby zobaczyć długość dyfuzji,
   średnie strumienie na ścianach oraz średni strumień \(J_b\) bez używania
   terminala. W `FrontAPP.py` liczba węzłów w każdym kierunku musi mieścić się
   w przedziale 3–50, wszystkie podawane wartości muszą być dodatnie, a
   porowatość \(\epsilon\) musi wynosić od 0,005 do 1. Współczynnik rozpadu
   \(\lambda\) jest ustawiony na stałą \(2{,}1\times10^{-6}\) 1/s, a w trakcie
   obliczeń pojawia się okno z paskiem postępu (aktualizowanym co 5%).
