#include <iostream>
#include <chrono>
#include <assert.h>
#include <getopt.h>

#include "model.hpp"
#include "benchmarks.cpp"

struct Domain {
  Domain(dotl::DoTL_model initial, 
         std::vector<dotl::Event_model> sequence_of_events, 
         std::unordered_map<dotl::Id, std::string> event_map) : 
    model(initial), sequence_of_events(sequence_of_events),
    event_map(event_map) {}
  dotl::DoTL_model model;
  std::vector<dotl::Event_model> sequence_of_events;
  std::unordered_map<dotl::Id, std::string> event_map;
};

// sets up domain for benchmarking:
// agents are {0,1,2}.
// the sequence of events is 5 flips, 1 mark, 7 flips, 1 lift.
// in the initial situation, the coin is unmarked and showing tails.
// flips: agent 0 and 1 takes turns secretly flipping (turning over) the coin.
// mark: agent 2 secretly marks the coin.
// lift: the marked coin showing tails is revealed to all agents.
Domain setup() {
  dotl::Formula T = dotl::Top();

  dotl::Event_model Flip0 = dotl::Event_model(3);
  dotl::Formula t = 't';
  dotl::Formula negt = dotl::Neg(t);
  auto flip0t = Flip0.new_event(true, negt, { 't' }, {});
  auto flip0h = Flip0.new_event(true, t, {}, { 't' });
  auto skip0 = Flip0.new_event(false, T, {}, {});
  Flip0.add_plausibility_relation(0, { {flip0t, flip0h} });
  Flip0.add_plausibility_relation(0, { {flip0h, flip0t} });
  Flip0.add_plausibility_relation(1, { {flip0t, skip0} });
  Flip0.add_plausibility_relation(1, { {flip0h, skip0} });
  Flip0.add_plausibility_relation(1, { {flip0t, flip0h} });
  Flip0.add_plausibility_relation(1, { {flip0h, flip0t} });
  Flip0.add_plausibility_relation(2, { {flip0t, skip0} });
  Flip0.add_plausibility_relation(2, { {flip0h, skip0} });
  Flip0.add_plausibility_relation(2, { {flip0t, flip0h} });
  Flip0.add_plausibility_relation(2, { {flip0h, flip0t} });

  dotl::Event_model Flip1 = dotl::Event_model(3);
  auto flip1t = Flip1.new_event(true, negt, { 't' }, {});
  auto flip1h = Flip1.new_event(true, t, {}, { 't' });
  auto skip1 = Flip1.new_event(false, T, {}, {});
  Flip1.add_plausibility_relation(0, { {flip1t, skip1} });
  Flip1.add_plausibility_relation(0, { {flip1h, skip1} });
  Flip1.add_plausibility_relation(0, { {flip1t, flip1h} });
  Flip1.add_plausibility_relation(0, { {flip1h, flip1t} });
  Flip1.add_plausibility_relation(1, { {flip1t, flip1h} });
  Flip1.add_plausibility_relation(1, { {flip1h, flip1t} });
  Flip1.add_plausibility_relation(2, { {flip1t, skip1} });
  Flip1.add_plausibility_relation(2, { {flip1h, skip1} });
  Flip1.add_plausibility_relation(2, { {flip1t, flip1h} });
  Flip1.add_plausibility_relation(2, { {flip1h, flip1t} });

  dotl::Event_model Mark = dotl::Event_model(3);
  dotl::Formula m = 'm';
  dotl::Formula negm = dotl::Neg(m);
  auto mark = Mark.new_event(true, negm, { 'm' }, {});
  auto skip = Mark.new_event(false, T, {}, {});
  Mark.add_plausibility_relation(0, { {mark, skip} });
  Mark.add_plausibility_relation(1, { {mark, skip} });

  dotl::Event_model ShowTailsAndMark = dotl::Event_model(3);
  dotl::Formula tm = dotl::Con(t, m);
  auto lift = ShowTailsAndMark.new_event(true, tm, {}, {});

  std::vector<dotl::Event_model> events;
  for (int i = 0; i < 12; i++) {
    if (i == 5) {
      events.emplace_back(Mark);
    }
    if (i % 2) {
      events.emplace_back(Flip1);
    }
    else {
      events.emplace_back(Flip0);
    };
  };
  events.emplace_back(ShowTailsAndMark);

  std::unordered_map<dotl::Id, std::string> event_map = 
      { 
            { flip0t, "flip0t" }, 
            { flip0h, "flip0h" },
            { skip0, "skip0" },
            { flip1t, "flip1t" }, 
            { flip1h, "flip1h" },
            { skip1, "skip1" },
            { mark, "mark" },
            { skip, "skip2" },
            { lift, "lift" }
      };

  dotl::DoTL_model D = dotl::DoTL_model(3, { 't' });

  return Domain(D, events, event_map);
};

struct Test_APU {
  Test_APU() : D(setup()) {
    next = D.sequence_of_events.begin();
  }
  void run() {
    if (next != D.sequence_of_events.end()) {
      D.model.action_priority_update(*next++);
    }
    //std::cerr << D.model.size() << " histories" << std::endl;
  };
  private:
  Domain D;
  std::vector<dotl::Event_model>::iterator next;
};

template <dotl::Abduction_type T>
struct Test_MPU {
  Test_MPU() : D(setup()) {
    for (auto& e: D.sequence_of_events) {
      e.perspective_shift(1);
    };
    next = D.sequence_of_events.begin();
  };
  void run() {
    if (degenerate)
      degenerate = !D.model.abduction_step(T, 1);
    else if (next != D.sequence_of_events.end()) {
      degenerate = !D.model.most_plausible_update(*next++, 1);
    }
    //std::cerr << D.model.size() << " histories" << std::endl;
  };
  private:
  Domain D;
  bool degenerate = false;
  std::vector<dotl::Event_model>::iterator next;
};

void make_benchmarks() {
  const int avg = 10;
  std::cerr << "generating benchmarks for most plausible update and "
            << "abduction with chronological minimization" << std::endl;
  benchmarks::Benchmark<Test_MPU<dotl::CHRON_MIN>, avg, 22> test_mpu_chron_min;
  test_mpu_chron_min.generate();
  std::cerr << test_mpu_chron_min.get_results();
  std::cerr << "generating benchmarks for most plausible update and "
            << "abduction with inverse chronological minimization" <<
            std::endl;
  benchmarks::Benchmark<Test_MPU<dotl::INV_CHRON_MIN>, avg, 22>
    test_mpu_inv_chron_min;
  test_mpu_inv_chron_min.generate();
  std::cerr << test_mpu_inv_chron_min.get_results();
  std::cerr << "generating benchmarks for most plausible update and "
            << "abduction with uniform expansion" << std::endl;
  benchmarks::Benchmark<Test_MPU<dotl::UNI_EXP>, avg, 17>
    test_mpu_uni_exp;
  test_mpu_uni_exp.generate();
  std::cerr << test_mpu_uni_exp.get_results();
  std::cerr << "generating benchmarks for action priority update" << std::endl;
  benchmarks::Benchmark<Test_APU, avg, 15> test_apu;
  test_apu.generate();
  std::cerr << test_apu.get_results();
};

void test_valuate(bool output_model = false) {
  Domain D = setup();
  std::unordered_map<std::string, dotl::Event_id> rev_event_map;
  for (auto & [k,v] : D.event_map) {
    rev_event_map[v] = k;
  };
  for (int i = 0; i < 2; i++) {
    D.model.action_priority_update(D.sequence_of_events[i]);
  };

  if (output_model) {
    D.model.to_dot(D.event_map, "flipflip_APU.dot");
  };
  
  dotl::Formula T = dotl::Top();
  dotl::Formula t = 't';
  dotl::Formula negt = dotl::Neg(t);

  dotl::Formula B2t = dotl::Believes(t, 2);
  dotl::Formula B1negt = dotl::Believes(negt, 1);
  dotl::Formula B0negt = dotl::Believes(negt, 0);
  dotl::Formula psi = dotl::Con(B1negt, B0negt);
  dotl::Formula phi = dotl::Con(B2t, psi);
  assert(D.model.valuate(phi, 3));
  //std::cout << "D, (w0;e" << rev_event_map["flip0h"] << ";e" 
  //          << rev_event_map["flip1t"] << ") |= " << phi 
  //          << " = " << D.model.valuate(phi, 3) << std::endl;

  dotl::Path path = 
    dotl::Seq(rev_event_map["flip0h"], rev_event_map["flip1t"]);
  dotl::Path inv_path = dotl::Inv(path);
  dotl::Formula right_history = dotl::After_event(T, inv_path);

  dotl::Path path1 = 
    dotl::Seq(rev_event_map["flip0h"], rev_event_map["skip1"]);
  dotl::Path inv_path1 = dotl::Inv(path1);
  dotl::Formula inner1 = dotl::After_event(T, inv_path1);
  dotl::Formula B0flipskip = dotl::Believes(inner1, 0);

  dotl::Path path2 = 
    dotl::Seq(rev_event_map["skip0"], rev_event_map["flip1h"]);
  dotl::Path inv_path2 = dotl::Inv(path2);
  dotl::Formula inner2 = dotl::After_event(T, inv_path2);
  dotl::Formula B1skipflip = dotl::Believes(inner2, 1);

  dotl::Path path3 = dotl::Seq(rev_event_map["skip0"], rev_event_map["skip1"]);
  dotl::Path inv_path3 = dotl::Inv(path3);
  dotl::Formula inner3 = dotl::After_event(T, inv_path3);
  dotl::Formula B2skipskip = dotl::Believes(inner3, 2);

  dotl::Formula con1 = dotl::Con(right_history, B0flipskip);
  dotl::Formula con2 = dotl::Con(B1skipflip, B2skipskip);
  dotl::Formula con3 = dotl::Con(con1, con2);
  assert (D.model.valuate(con3, 3));
  //std::cout << "D, (w0;e" << rev_event_map["flip0h"] << ";e" 
  //          << rev_event_map["flip1t"] << ") |= " << con3 << " = " 
  //          << D.model.valuate(con3, 3) << std::endl;
};

void test_abduction(bool output_models = false) {
  Domain D1 = setup();

  for (int i = 0; i < 2; i++) {
    D1.model.action_priority_update(D1.sequence_of_events[i]);
  };
  dotl::Event_model ShowTails = dotl::Event_model(3);
  dotl::Formula t = 't';
  auto lift = ShowTails.new_event(true, t, {}, {});
  D1.model.action_priority_update(ShowTails);
  D1.event_map[lift] = "lift";
  if (output_models) {
    D1.model.to_dot(D1.event_map, "flipfliplift_APU.dot");
  };

  Domain D2 = setup();
  for (int i = 0; i < 2; i++) {
    D2.sequence_of_events[i].perspective_shift(1);
    assert(D2.model.most_plausible_update(D2.sequence_of_events[i], 1));
  };
  assert(
    !D2.model.most_plausible_update(ShowTails, 1));
  if (output_models) {
    D2.model.to_dot(D2.event_map, "flipfliplift_MPU.dot");
  };

  assert(D2.model.abduction_step(dotl::CHRON_MIN, 1));
  if (output_models) {
    D2.model.to_dot(D2.event_map, "flipfliplift_ABDUCTED.dot");
  };
};

int main(int argc, char* argv[]) {
  int opt;
  while ((opt = getopt(argc, argv, "tbh")) != -1) {
    switch (opt) {
      case 't':
        std::cerr << "testing abduction..." << std::endl;
        test_abduction(true);
        std::cerr << "testing valuate... " << std::endl;
        test_valuate(true);
        std::cerr << "Everything ok. Models can be compiled to pdf by " <<
        "running $ dot -Tps <name.dot> -o <name.pdf>" << std::endl;
        break;
      case 'b':
        make_benchmarks();
        break;
      case 'h':
        std::cerr << "Usage: " << argv[0] << " [-tbh]" << std::endl;
        std::cerr << "  -t           tests abduction and valuate, generating" << std::endl;
        std::cerr << "               .dot files for the models used. These" << std::endl;
        std::cerr << "               can be compiled to pdf using:" << std::endl;
        std::cerr << "               $ dot -Tps <name.dot> -o <name.pdf>" << std::endl;
        std::cerr << "  -b           generates benchmarks." << std::endl;
        std::cerr << "  -h           displays this message." << std::endl;
        break;
    };
  };
  if (argc == 1) {
    std::cerr << "Usage: " << argv[0] << " [-tbh]" << std::endl;
    std::cerr << "  -t           tests abduction and valuate, generating" << std::endl;
    std::cerr << "               .dot files for the models used. These" << std::endl;
    std::cerr << "               can be compiled to pdf using:" << std::endl;
    std::cerr << "               $ dot -Tps <name.dot> -o <name.pdf>" << std::endl;
    std::cerr << "  -b           generates benchmarks." << std::endl;
    std::cerr << "  -h           displays this message." << std::endl;
    return 0;
  };
}
